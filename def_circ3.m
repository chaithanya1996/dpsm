function [pos,norms]=def_circ3(dr, rad, center, normal)
%def_circ3 円盤状音源の定義
% [pos,norms]=def_circ3(dr, rad, center, normal)
% 引数
% dr : 音源間の距離
% rad: 円盤半径
% center: 円盤の中心
% normal: 円盤の法線 (放射方向)
%
% 戻り値 [音源数, 次元数=3]の行列
% pos:  音源位置
% norm: 音源の法線ベクトル = normalの音源数コピー

n = normal/norm(normal);
l = cross(normal, [0 1 0]);
l = l/norm(l);
m = cross(n, l);
rot = [l;m;n];

pos = []; norms=[];
for ri=dr:dr:rad
    an = ceil(2*pi*ri/dr);
    for ai=0:an-1
        dpos = ri* [cos(2*pi*ai/an); sin(2*pi*ai/an); 0];
        npos = center' + rot'*dpos;
        pos = [pos; npos'];
        norms = [norms; normal];
    end
end

end
