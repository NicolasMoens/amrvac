;.r getmflow
;askstr,'func(s) (e.g. rho p m1 Te;m2 r*m1 -T) ',func,doask
func='m1'
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    ;read,cp,prompt='Enter the position of center point: '
    cpl=6.
    cpr=20.
    ;if(cp lt domainamr[0] or cp gt domainamr[1]) then message,'The center point is out of domain!'
    ;read,radius,prompt='Enter radius of the range: '
    radius=0.
    ;radius=abs(radius)
    if(radius eq 0.) then begin
      av1=getvalue(cpl,func,xamr,wamr,ngrids,ndim,nx,nw,physics,eqpar,wnames)
      av2=getvalue(cpr,func,xamr,wamr,ngrids,ndim,nx,nw,physics,eqpar,wnames)
      av=av1-av2
    endif else begin
      if(cp-radius lt domainamr[0] or cp+radius gt domainamr[1]) then message,'The radius is too large!'
      read,npt,prompt='Enter the number of points in the range: '
      dx=2.*radius/float(npt-1)
      av=0.
      for i=0,npt-1 do begin
        xp=cp-radius+i*dx
        av=av+getvalue(xp,func,xamr,wamr,ngrids,ndim,nx,nw,physics,eqpar,wnames)
      endfor
    av=av/float(npt)
    endelse
  end
  ;-------------- 2D ----------------------
  2: begin
    read,cp1,prompt='Enter the position of center point x: '
    read,cp2,prompt='Enter the position of center point y: '
    if(cp1 lt domainamr[0] or cp1 gt domainamr[2]) then message,'x of center point is out of domain!'
    if(cp2 lt domainamr[1] or cp2 gt domainamr[3]) then message,'y of center point is out of domain!'
    read,radius,prompt='Enter radius of the range: '
    xp=fltarr(ndim)
    radius=abs(radius)
    if(radius eq 0.) then begin
      xp[0]=cp1
      xp[1]=cp2
      av=getvalue(xp,func,xamr,wamr,ngrids,ndim,nx,nw,physics,eqpar,wnames)
    endif else begin
      if(cp1-radius lt domainamr[0] or cp1+radius gt domainamr[2]) then message,'The radius is too large for x direction!'
      if(cp2-radius lt domainamr[1] or cp2+radius gt domainamr[3]) then message,'The radius is too large for y direction!'
      read,npt,prompt='Enter the number of points in the range: '
      dx=2.*radius/float(npt-1)
      av=0.
      for j=0,npt-1 do begin
        for i=0,npt-1 do begin
          xp[0]=cp1-radius+i*dx
          xp[1]=cp2-radius+j*dx
          av=av+getvalue(xp,func,xamr,wamr,ngrids,ndim,nx,nw,physics,eqpar,wnames)
        endfor
      endfor
      av=av/float(npt*npt)
    endelse
  end
  ;-------------- 3D ----------------------
  3: begin
    read,cp1,prompt='Enter the position of center point x: '
    read,cp2,prompt='Enter the position of center point y: '
    read,cp3,prompt='Enter the position of center point z: '
    if(cp1 lt domainamr[0] or cp1 gt domainamr[3]) then message,'x of center point is out of domain!'
    if(cp2 lt domainamr[1] or cp2 gt domainamr[4]) then message,'y of center point is out of domain!'
    if(cp3 lt domainamr[2] or cp2 gt domainamr[5]) then message,'z of center point is out of domain!'
    read,radius,prompt='Enter radius of the range: '
    xp=fltarr(ndim)
    radius=abs(radius)
    if(radius eq 0.) then begin
      xp[0]=cp1
      xp[1]=cp2
      xp[2]=cp3
      av=getvalue(xp,func,xamr,wamr,ngrids,ndim,nx,nw,physics,eqpar,wnames)
    endif else begin
      if(cp1-radius lt domainamr[0] or cp1+radius gt domainamr[3]) then message,'The radius is too large for x direction!'
      if(cp2-radius lt domainamr[1] or cp2+radius gt domainamr[4]) then message,'The radius is too large for y direction!'
      if(cp3-radius lt domainamr[2] or cp2+radius gt domainamr[5]) then message,'The radius is too large for z direction!'
      read,npt,prompt='Enter the number of points in the range: '
      dx=2.*radius/float(npt-1)
      av=0.
      xp=fltarr(ndim)
      for k=0,npt-1 do begin
      for j=0,npt-1 do begin
      for i=0,npt-1 do begin
        xp[0]=cp1-radius+i*dx
        xp[1]=cp2-radius+j*dx
        xp[2]=cp3-radius+k*dx
        av=av+getvalue(xp,func,xamr,wamr,ngrids,ndim,nx,nw,physics,eqpar,wnames)
      endfor
      endfor
      endfor
      av=av/float(npt*npt*npt)
    endelse
  end
  endcase
  print,'average value of '+'func: ',av
end
