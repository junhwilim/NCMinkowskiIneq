# NCMinkowskiIneq
MATLAB codes for computations in "J. Lim, Noncommutative Minkowski integral inequality and a unitary categorification criterion for fusion rings, 2026. arXiv:2601.13490."<br>
To load computation result, run `load('FusionRingDataSet_Solutions.mat', 'FusionRingDataSet_Solutions')`. <br>
* `FusionRingDataSet_Solutions{m}{r}{n}`: `n`th fusion ring with multiplicity `m` and rank `r+1`<br>
* `FusionRingDataSet_Solutions{m}{r}{n}{l}`: `l`th solution of RHS/LHS<1 <br>
* `FusionRingDataSet_Solutions{m}{r}{n}{l}{t}` with `t`=1,2,3: the index of basis element <br>
* `FusionRingDataSet_Solutions{m}{r}{n}{l}{4}`: $$((a_n)_{n=1}^{r+1}, 1/p)$$ <br>
* `FusionRingDataSet_Solutions{m}{r}{n}{l}{5}`: RHS/LHS <br>

To run a new computation, run `FusionRingDataSet0.m`.
