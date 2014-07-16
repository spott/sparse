{-# LANGUAGE MultiParamTypeClasses #-}

module Data.Sparse where

import qualified Data.Vector.Unboxed as U

type Index = Integer


class (Num val) => Vector v val where
        (<!!>) :: v val -> Int -> val

class (Num val) => Matrix m val where
        -- multiply a matrix by a vector
        mul :: (Vector v val) => m val  -> v val -> v val
        -- update the diagonal
        updateDiagonal :: m val -> (Int -> val) -> m val
        -- get a row:
        (<!>) :: m val -> Int -> U.Vector val

data CSparseRow val = CSparseRow     { values     :: U.Vector val
                                     , colIndexes :: U.Vector Int
                                     --the start index in values and
                                     --colIndices of each row.  should be
                                     --terminated with "rows"
                                     , rowIndexes :: U.Vector Int
                                     , rows       :: Int
                                     , cols       :: Int }
                                     deriving Show

csparseRowGetRow :: (Num val) => CSparseRow val -> Int -> U.Vector val
csparseRowGetRow m row = let vec = U.replicate c 0 in
                             U.update_ vec (part (colIndexes m)) (part (values m))
        where part = U.slice (rowIndexes m U.! row) (rowIndexes m U.! (row + 1))
              r = rows m
              c = cols m

{-CSparseRowMul :: (Vector v ord val) => CSparseRow ord val -> v ord val -> v ord val-}
{-CSparseRowMul m v = -}

--instance Matrix m ord val => Sparse m ord val where


