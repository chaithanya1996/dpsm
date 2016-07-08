function [pos, norms]=def_line3(dr, len, center, normal)
%def_line3 線状音源の定義
% [pos,norms]=def_line3(dr, len, center, normal)
% 引数
% dr : 音源間の距離
% len: 線音源長
% center: 線音源の中心
% normal: 線音源の法線 (放射方向)
%
% 戻り値 [音源数, 次元数=3]の行列
% pos:  音源位置
% norm: 音源の法線ベクトル = normalの音源数コピー

n = normal/norm(normal);
l = cross(normal, [0 0 1]);
l = l/norm(l);
m = cross(n, l);
rot = [l;m;n];


XI=[-len/2:dr:len/2];

pnum = length(XI(:));
dpos = [XI(:), zeros(pnum,1), zeros(pnum,1)];
rpos = rot'*dpos';

pos = (rpos + center'*ones(1,pnum))';
norms = ones(pnum,1)*normal;

end
