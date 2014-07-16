{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Data.Sparse where

import qualified Data.Vector.Unboxed as U

type Index = Integer


{-class (Num val) => Vector v val where-}
        {-(!) :: v val -> Int -> val-}

class (Num val) => Matrix m val where
        -- multiply a matrix by a vector
        (<**>) :: m val -> U.Vector val -> U.Vector val
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
                                     , rows       :: !Int
                                     , cols       :: !Int }
                                     deriving Show

{-diagonalMatrix :: (U.Unbox val, Num val) => (Int, Int) -> (Int -> val) -> CSparseRow val-}
{-diagonalMatrix size f = -}

csparseRowGetRow :: (U.Unbox val, Num val) => CSparseRow val -> Int -> U.Vector val
csparseRowGetRow m row = let vec = U.replicate c 0 in
                             U.update_ vec partcolIndexes partValues
        where partValues = U.slice (rowIndexes m U.! row) (rowIndexes m U.! (row + 1)) (values m)
              partcolIndexes = U.slice (rowIndexes m U.! row) (rowIndexes m U.! (row + 1)) (colIndexes m)
              c = cols m

csparseRowMul :: (U.Unbox val, Num val) => CSparseRow val -> U.Vector val -> U.Vector val
csparseRowMul m v = U.generate r rowCalc
        where rowCalc row = U.sum $ U.ifilter (\ i _ -> i `U.elem` partcolIndexes row) $ U.accumulate_ (*) v (partcolIndexes row) (partValues row)
              partValues row = U.slice (slice row) (slice (row + 1) - slice row) (values m)
              partcolIndexes row = U.slice (slice row) (slice (row + 1) - slice row) (colIndexes m)
              slice row = U.snoc (rowIndexes m) r U.! row
              r = rows m

csparseUpdateDiagonal :: (Num val) => CSparseRow val -> (Int -> val) -> CSparseRow val
csparseUpdateDiagonal m f = m

instance (U.Unbox val, Num val) => Matrix CSparseRow val where
        (<**>) = csparseRowMul
        (!!) = csparseRowGetRow
        updateDiagonal = csparseUpdateDiagonal

{-instance (Num val) => Vector U.Vector val where-}
        {-(!) = (U.!)-}


--
