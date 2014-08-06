{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Data.Sparse where

import qualified Data.List           as L
import           Debug.Trace
--import qualified Data.Vector.Generic as V
import qualified Data.Vector.Unboxed as U
import           Prelude             hiding (pred, (!!))

type Index = (Int, Int)
type Size = Index

class (Num val, Pattern p) => MatrixClass m p val where
        -- multiply a matrix by a vector
        (*:) :: m p val -> U.Vector val -> U.Vector val
        generateMatrix :: (Num val, Pattern p) => p -> (Index -> val) -> m p val
        -- get a row:
        (!!) :: m p val -> Int -> U.Vector val
        {-(+) :: m p val -> m p val -> m p val-}
        {-(-) :: m p val -> m p val -> m p val-}

class Pattern p where
        union :: p -> p -> p
        eye :: Size -> p
        index :: p -> Int -> Index
        size :: p -> Size
        generatePattern :: Size -> (Index -> Bool) -> p
        nnz :: p -> Int

data CompressedSparseRow = CompressedSparseRow { _colIndexes   :: !(U.Vector Int)
                                                 --the start index in values and
                                                 --colIndices of each row.
                                                 , _rowIndexes :: !(U.Vector Int)
                                                 --if length rowIndexes
                                                 -- is < (fst size), then
                                                 -- there are empty rows.
                                                 , _sizeCSR    :: !Size} deriving Show

data RowMajorDense = RowMajorDense { _sizeRMD :: !Index } deriving Show

instance Pattern RowMajorDense where
        eye s = RowMajorDense { _sizeRMD = s }
        index p = \i -> (i `div` snd (size p), i `mod` snd (size p))
        union p1 p2 | size p1 == size p2 = RowMajorDense { _sizeRMD = size p1 }
                    | otherwise = undefined
        size = _sizeRMD
        generatePattern s _ = RowMajorDense { _sizeRMD = s }
        nnz p = (fst $ size p) * (snd $ size p)

instance Pattern CompressedSparseRow where
        eye s = CompressedSparseRow { _colIndexes = U.generate (min (fst s) (snd s)) (\i -> i)
                                    , _rowIndexes = U.generate (min (fst s) (snd s)) (\i -> i)
                                    , _sizeCSR = s }

        index p = \i -> let c = _colIndexes p U.! i
                            rs = _rowIndexes p
                       in  (U.ifoldl' ( \le i' e -> if e <= i then i' else le) (fst $ size p) rs, c)

        union p1 p2 | size p1 /= size p2 = undefined
                    | size p1 == size p2 = CompressedSparseRow { _colIndexes = U.concat colslist
                                          , _rowIndexes =  U.scanl' (+) 0 $ U.fromList $ map U.length colslist
                                          , _sizeCSR = size p1}
                      where
                            colslist = map (\ row -> mergeUnion (partcolIndexes row p1) (partcolIndexes row p2)) [0..(r p1)]
                            partcolIndexes row p = U.slice (slice p row) (slice p (row + 1) - slice p row) (_colIndexes p)
                            slice p = (U.!) (U.snoc (_rowIndexes p) $ r p)
                            r p = fst $ size p

        size = _sizeCSR

        generatePattern s pred = CompressedSparseRow { _colIndexes = L.foldl1' (U.++) filteredvalues
                                              , _rowIndexes = U.prescanl (+) 0 $  U.fromList $ map U.length $ filteredvalues
                                              , _sizeCSR = s}
                        where filteredvalues = map (\i -> snd $ U.unzip (U.filter pred i)) $ allvalues
                              allvalues = map (\r -> U.zip (U.replicate cols r) (U.enumFromN 0 cols)) $ [0..rows-1]
                              rows = fst s
                              cols = snd s

        nnz = U.length . _colIndexes

traced :: (Show a) => a -> a
traced a = traceShow a a

mergeUnion :: (U.Unbox val, Ord val) => U.Vector val -> U.Vector val -> U.Vector val
mergeUnion a b | U.null a = b
               | U.null b = a
               | (U.head a) == (U.head b) = (U.head a) `U.cons` mergeUnion (U.tail a) (U.tail b)
               | (U.head a) < (U.head b) = (U.head a) `U.cons` mergeUnion (U.tail a) b
               | (U.head a) > (U.head b) = (U.head b) `U.cons` mergeUnion a (U.tail b)
               | otherwise = undefined

data Matrix p val = Matrix { _values :: !(U.Vector val), _pattern :: p }

instance (Show val, U.Unbox val, Num val) => MatrixClass Matrix CompressedSparseRow val where
        (*:) m v = U.generate r rowCalc
            where rowCalc row = U.sum $ U.ifilter (\ i _ -> i `U.elem` partcolIndexes row) $ U.accumulate_ (*) v (partcolIndexes row) (partValues row)
                  partValues row = U.slice (slice row) (slice (row + 1) - slice row) (_values m)
                  partcolIndexes row = U.slice (slice row) (slice (row + 1) - slice row) (_colIndexes pat)
                  slice row = U.snoc (_rowIndexes pat) (nnz pat) U.! row
                  r = fst $ size pat
                  pat = _pattern m
        generateMatrix p f = Matrix { _values = U.generate (nnz p) $ f . (index p), _pattern = p }

        (!!) m row = let vec = U.replicate c 0 in
                                 U.update_ vec partcolIndexes partValues
                where
                      partValues = U.slice (slice row) (slice (row + 1) - slice row) (_values m)
                      partcolIndexes = U.slice (slice row) (slice (row + 1) - slice row) (_colIndexes p)
                      slice = (U.!) (U.snoc (_rowIndexes p) $ nnz p)
                      c = snd $ size p
                      p = _pattern m
        {-(+) m1 m2 = Matrix { _values = -}
                           {-, _pattern = newpattern}-}
                           {-where v1 = generateMatrix newpattern (\(r,c) -> )-}
                                 {-newpattern = union (_pattern m1) (_pattern m2)-}

instance (Show val, U.Unbox val, Num val) => Show (Matrix CompressedSparseRow val) where
        show m = concatMap (\x -> showVec (m !! x) ++ "\n") [0..fst (size $ _pattern m) - 1]
            where showVec v = U.foldl' (\a b -> a ++ ", " ++ show b) ("[" ++ show (U.head v)) (U.tail v) ++ "]"

{---This might want to be in a function...-}
{-imap :: (U.Unbox a, U.Unbox b) => (Index -> a -> b) -> CSparseRow a -> CSparseRow b-}
{-imap f m = m {values =  U.map (uncurry f) (U.zip index (values m)) }-}
        {-where index :: U.Vector Index-}
              {-index = U.zip rows (colIndexes m)-}
              {-rows = U.convert $ V.foldl1' (V.++) $ V.zipWith (\a b -> V.replicate (b - a) a) ri (V.drop 1 ri)-}
              {-ri = V.convert $ U.snoc (rowIndexes m) (fst $ size m)-}
