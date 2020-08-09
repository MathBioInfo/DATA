function stampname(xpos,ypos,rot)

dn = ['[',date,' Amjad Khan, Alita R. Burmeister, Lindi M. Wahl, Evolution along the parasitism-mutualism continuum determines the genetic repertoire of prophages]'];
tpos = text(xpos,ypos,dn,'fontsize',10);
set(tpos,'rotation',rot);
