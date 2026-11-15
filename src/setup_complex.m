function cmplx = setup_complex(dim)
% Empty cell complex of dimension dim.
cmplx(dim).name='';

cmplx(dim).num(dim).val=[];
cmplx(dim).break=[];

cmplx(dim).cc=[];
cmplx(dim).bc=[];

cmplx(dim).dir=[];

cmplx(dim).vvol=[];
cmplx(dim).dvol=[];
cmplx(dim).bfvol=[];
cmplx(dim).fvol=[];

cmplx(dim).omega=[];

cmplx(dim).bndop(dim).indx=[];
cmplx(dim).bndop(dim).sgn=[];

cmplx(dim).flat=[];
cmplx(dim).sharp=[];
