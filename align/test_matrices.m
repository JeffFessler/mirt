ig = image_geom('nx', 50, 'ny', 50, 'dx', 1, 'dy', 1);
kg = knot_geom('nx', 10, 'ny', 10, 'mx', 5, 'my', 5);

clear Y Y2;

[B Bgx Bgy] = makeB(ig, kg);

I = eye(100);

for i = 1 : 100
	Y(:, i) = B * I(:, i);
end

I2 = eye(2500);

for i = 1 : 2500
	Y2(:, i) = B' * I2(:, i);
end

printf('Accuracy of differences for B = %d\n'...
	, max(max(abs(Y2-Y'))) / max(max(abs(Y2))));

for i = 1 : 100
	Y(:, i) = Bgx * I(:, i);
end

for i = 1 : 2500
	Y2(:, i) = Bgx' * I2(:, i);
end

printf('Accuracy of differences for Bgx = %d\n'...
	, max(max(abs(Y2-Y'))) / max(max(abs(Y2))));


for i = 1 : 100
	Y(:, i) = Bgy * I(:, i);
end

for i = 1 : 2500
	Y2(:, i) = Bgy' * I2(:, i);
end

printf('Accuracy of differences for Bgy = %d\n'...
	, max(max(abs(Y2-Y'))) / max(max(abs(Y2))));


for i = 1 : 2500
	Y2(:, i) = Bgy' * I2(:, i);
end

printf('Accuracy of differences for Bgx/Bgy = %d\n'...
	, max(max(abs(Y2-Y'))) / max(max(abs(Y2))));

A1 = single(rand([10 10]));
A2 = single(rand([10 10]));



[Bw Bwgx Bwgy] = makeB(ig, kg, B*A1, B*A2);

for i = 1 : 100
	Y(:, i) = Bw * I(:, i);
end

I2 = eye(2500);

for i = 1 : 2500
	Y2(:, i) = Bw' * I2(:, i);
end

printf('Accuracy of differences for Bw = %d\n'...
	, max(max(abs(Y2-Y'))) / max(max(abs(Y2))));


for i = 1 : 100
	Y(:, i) = Bwgx * I(:, i);
end

for i = 1 : 2500
	Y2(:, i) = Bwgx' * I2(:, i);
end

printf('Accuracy of differences for Bwgx = %d\n' ...
	, max(max(abs(Y2-Y'))) / max(max(abs(Y2))));



for i = 1 : 2500
	Y2(:, i) = Bwgy' * I2(:, i);
end

printf('Accuracy of differences for Bwgx/Bwgy = %d\n'...
	, max(max(abs(Y2-Y'))) / max(max(abs(Y2))));




[W Wgx Wgy] = makeW({B, B}, {A1, A2});

clear Y Y2;

for i = 1 : 2500
	Y(:, i) = W * I2(:, i);
end

for i = 1 : 2500
	Y2(:, i) = W' * I2(:, i);
end

printf('Accuracy of differences for W = %d\n'...
	, max(max(abs(Y2-Y'))) / max(max(abs(Y2))));


for i = 1 : 2500
	Y(:, i) = Wgx * I2(:, i);
end

for i = 1 : 2500
	Y2(:, i) = Wgx' * I2(:, i);
end

printf('Accuracy of differences for Wgx = %d\n'...
	, max(max(abs(Y2-Y'))) / max(max(abs(Y2))));


for i = 1 : 2500
	Y2(:, i) = Wgy' * I2(:, i);
end

printf('Accuracy of differences for Wgx/Wgy = %d\n'...
	, max(max(abs(Y2-Y'))) / max(max(abs(Y2))));

