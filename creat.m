function creat(name)
nameps = sprintf('%s.eps',name);
print(nameps,'-depsc2');
namjpg = sprintf('%s.jpg',name);
print(namjpg,'-djpeg80');

