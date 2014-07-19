{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Data.Sparse where

import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U
import           Prelude             hiding ((!!))

type Index = (Int, Int)

{-class (Num val) => Vector v val where-}
        {-(!) :: v val -> Int -> val-}

class (Num val) => Matrix m val where
        -- multiply a matrix by a vector
        (**) :: m val -> U.Vector val -> U.Vector val
        generate :: (U.Unbox val, Num val) => Index -> (Index -> Bool) -> (Index -> val) -> m val
        -- update the diagonal
        updateDiagonal :: m val -> (Int -> val) -> m val
        -- get a row:
        (!!) :: m val -> Int -> U.Vector val
        dim :: m val -> Index

data CSparseRow val = CSparseRow     { values     :: !(U.Vector val)
                                     , colIndexes :: !(U.Vector Int)
                                     --the start index in values and
                                     --colIndices of each row.  should be
                                     --terminated with "rows"
                                     , rowIndexes :: !(U.Vector Int)
                                     , size       :: !Index}

instance (Show val, U.Unbox val, Num val) => Show (CSparseRow val) where
        show = showMatrix
--Slow
showMatrix :: (Show val, U.Unbox val, Num val) => CSparseRow val -> String
showMatrix m = concatMap (\x -> showVec (m !! x) ++ "\n") [0..fst (dim m) - 1]
        where showVec v = U.foldl' (\a b -> a ++ ", " ++ show b) ("[" ++ show (U.head v)) (U.tail v) ++ "]"

diagonalMatrix :: (Matrix m val, U.Unbox val, Num val) => Index -> (Int -> val) -> m val
diagonalMatrix s f = generate s (uncurry (==)) (\(x,_) -> f x)

generateCSparseRow :: (U.Unbox val, Num val) => Index -> (Index -> Bool) -> (Index -> val) -> CSparseRow val
generateCSparseRow s p func = CSparseRow { values = U.fromList $ foldl1 (++) (map vals [0..(fst s)]),
                                       colIndexes = U.fromList $ foldl1 (++) (map cols [0..(fst s)]),
                                       rowIndexes = U.fromList rowindex,
                                       size = s}
        where cols row = [col | col <- [0..(snd s)], (\i -> p (row,i)) col]
              vals row = [(\i -> func (row,i)) col | col <- cols row]
              rowindex = scanl (\row l -> l + length (cols row)) 0 [0..(fst s)]

rowOfCSparseRow :: (U.Unbox val, Num val) => CSparseRow val -> Int -> U.Vector val
rowOfCSparseRow m row = let vec = U.replicate c 0 in
                             U.update_ vec partcolIndexes partValues
        where
              partValues = U.slice (slice row) (slice (row + 1) - slice row) (values m)
              partcolIndexes = U.slice (slice row) (slice (row + 1) - slice row) (colIndexes m)
              slice = (U.!) (U.snoc (rowIndexes m) r)
              c = snd $ size m
              r = fst $ size m

mulCSparseRow :: (U.Unbox val, Num val) => CSparseRow val -> U.Vector val -> U.Vector val
mulCSparseRow m v = U.generate r rowCalc
        where rowCalc row = U.sum $ U.ifilter (\ i _ -> i `U.elem` partcolIndexes row) $ U.accumulate_ (*) v (partcolIndexes row) (partValues row)
              partValues row = U.slice (slice row) (slice (row + 1) - slice row) (values m)
              partcolIndexes row = U.slice (slice row) (slice (row + 1) - slice row) (colIndexes m)
              slice row = U.snoc (rowIndexes m) r U.! row
              r = fst $ size m

--This might want to be in a function...
imap :: (U.Unbox a, U.Unbox b) => (Index -> a -> b) -> CSparseRow a -> CSparseRow b
imap f m = m {values =  U.map (uncurry f) (U.zip index (values m)) }
        where index :: U.Vector Index
              index = U.zip rows (colIndexes m)
              rows = U.convert $ V.foldl1' (V.++) $ V.zipWith (\a b -> V.replicate (b - a) a) ri (V.drop 1 ri)
              ri = V.convert $ U.snoc (rowIndexes m) (fst $ size m)

csparseUpdateDiagonal :: (Num val) => CSparseRow val -> (Int -> val) -> CSparseRow val
csparseUpdateDiagonal m f = m

instance (U.Unbox val, Num val) => Matrix CSparseRow val where
        (**) = mulCSparseRow
        (!!) = rowOfCSparseRow
        generate = generateCSparseRow
        updateDiagonal = csparseUpdateDiagonal
        dim = size


--
