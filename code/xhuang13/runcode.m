%testing for rmin
for i=3:2:10
  mkdir(['rmin1/', num2str(i)]);
  folder = ['rmin1/', num2str(i), '/'];
  mainFx(20, 40, 0.4, i, 1, folder, 40, 40);
end

%testing for effect macro dimension

%testing for 