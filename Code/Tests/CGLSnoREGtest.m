%%Quick test to make sure CGLS_No REg Solves the equation

A = [1,24,3;4,10,6;17,28,9];

d = [9;3;4];

OPERATOR = @Optest;
PARAM.A = A;

[x] = cgls_o_noReg(OPERATOR, PARAM, [0;0;0], d, 1000, 0.005);
