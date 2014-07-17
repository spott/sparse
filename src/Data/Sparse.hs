{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Data.Sparse where

import qualified Data.Vector.Unboxed as U

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



data CSparseRow val = CSparseRow     { values     :: !(U.Vector val)
                                     , colIndexes :: !(U.Vector Int)
                                     --the start index in values and
                                     --colIndices of each row.  should be
                                     --terminated with "rows"
                                     , rowIndexes :: !(U.Vector Int)
                                     , size :: !Index}
                                     deriving Show

diagonalMatrix :: (Matrix m val, U.Unbox val, Num val) => Index -> (Int -> val) -> m val
diagonalMatrix size f = generate size (\(x,y) -> x == y) (\(x,y) -> f x)

generateCSparseRow :: (U.Unbox val, Num val) => Index -> (Index -> Bool) -> (Index -> val) -> CSparseRow val
generateCSparseRow size pred func = CSparseRow { values = U.fromList $ foldl1 (++) (map vals [0..(fst size)]), 
                                       colIndexes = U.fromList $ foldl1 (++) (map cols [0..(fst size)]), 
                                       rowIndexes = U.fromList $ rowindex,
                                       size = size}
        where cols row = [i | i <- [0..(snd size)], (\i -> pred $ (row,i)) i]
              vals row = [(\i -> func $ (row,i)) i | i <- (cols row)]
              rowindex = scanl (\row last -> last + (length $ cols row)) 0 [0..(fst size)]


rowOfCSparseRow :: (U.Unbox val, Num val) => CSparseRow val -> Int -> U.Vector val
rowOfCSparseRow m row = let vec = U.replicate c 0 in
                             U.update_ vec partcolIndexes partValues
        where partValues = U.slice (rowIndexes m U.! row) (rowIndexes m U.! (row + 1)) (values m)
              partcolIndexes = U.slice (rowIndexes m U.! row) (rowIndexes m U.! (row + 1)) (colIndexes m)
              c = snd $ size m

mulCSparseRow :: (U.Unbox val, Num val) => CSparseRow val -> U.Vector val -> U.Vector val
mulCSparseRow m v = U.generate r rowCalc
        where rowCalc row = U.sum $ U.ifilter (\ i _ -> i `U.elem` partcolIndexes row) $ U.accumulate_ (*) v (partcolIndexes row) (partValues row)
              partValues row = U.slice (slice row) (slice (row + 1) - slice row) (values m)
              partcolIndexes row = U.slice (slice row) (slice (row + 1) - slice row) (colIndexes m)
              slice row = U.snoc (rowIndexes m) r U.! row
              r = fst $ size m

csparseUpdateDiagonal :: (Num val) => CSparseRow val -> (Int -> val) -> CSparseRow val
csparseUpdateDiagonal m f = m

instance (U.Unbox val, Num val) => Matrix CSparseRow val where
        (**) = mulCSparseRow
        (!!) = rowOfCSparseRow
        generate = generateCSparseRow
        updateDiagonal = csparseUpdateDiagonal

{-instance (Num val) => Vector U.Vector val where-}
        {-(!) = (U.!)-}


--
