function data = get_ifndata(r, rows)
if r == 14
    nt =[[0.15,    0.3,    0.6,    1.2,    2.375,  3,      4.75]; ...               % nt yields dose of NP
              [0.94,    1.875,  3.75,   7.5,    NaN,    NaN,    NaN];...
              [0.19,    0.375,  0.75,   1.5,    3,      6,      NaN];...
              [1.5,     3,      6,      12,     24,     48,     NaN];...
              [2.5,     5,      10,     20,     40,     NaN,    NaN]];

    ifn = [[0.09,  0.19,   0.29,   0.41,   0.976,  1.9,    2.27]; ...               % Actual ifng data points
                [0.06,  0.13,   0.49,   0.85,   NaN,    NaN,    NaN];...
                [0.015, 0.021,  0.025,  0.055,  0.143,  0.241,  NaN];...
                [0.2,   0.31,   0.5,    1.09,   2.32,   2.9,    NaN];...
                [0.01,  0.03,   0.1,    0.29,    0.75,  NaN,    NaN]];
     v = [54; 31; 14; 11; 8];                                                       % 5 'rows' for r=14 corresponding to each valence
elseif r==20
    nt =[[2.0,   1.0,     0.5,      0.25,    0.125,   0.061]; ...
              [4.0,     2.0,    1.0,    0.5,    0.25,    0.125];...
              [75,      37.5,  	19,     10,     5,      NaN];...
              [4.0,     2.0,    1.0,    0.5,      NaN,    NaN]];

    ifn = [[1.86,   0.88,   0.26,   0.12,   0.07,   0.04]; ...
                [1.05,	0.41,	0.2,	0.085,	0.03,	0.01];...
                [0.279,	0.166,	0.07,	0.03,	0.004,  NaN];...
                [0.04	0.02	0.005	0.001,  NaN,    NaN]];
     v = [210; 61; 13; 9];
end
if nargin ==1
    rows = 1:size(ifn,1);
end
data = [];
for i = rows
  for j = 1:size(ifn,2)
    if isnan(ifn(i,j)) == 0
        data = [data;[v(i),nt(i,j),ifn(i,j)]];
    end
  end
end
end

