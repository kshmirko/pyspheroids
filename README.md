# pyspheroids


# Состав пакета
This module 'libspheroids' is auto-generated with f2py (version:2).

### Functions:


  alloc_dls_array(key,keyel,key_alloc)
  
	matrix_fix(kn1,grid1,wavel,kre,kim,are,aim,ratio,ndp)
  
	dlnr1 = usmatrix(ndp)
  
	usuf_ls(kn1,grid1,wavel,kre,kim,are,aim)
  
	usu_ls(ratio,nratn,kn1,grid1,wavel,kre,kim,are,aim)
  
	usu_ls_rd(ratio,nratn,kn1,grid1,wavel,kre,kim,are,aim)
  
	optchar(ndp)
  
	dls_read_input(fullpath_input)
  
	rrr,ar,ac = sizedisdn(kn,ia,id,cm,sm,rmm,rmin,rmax,knpar,idbg,nmd=len(cm))
  
	ac = sdnorm(c,s,rm,rmin,rmax)
  
	slog = slog(c,s,rm,rr)
  
	sizedisdn1(kn,ia,id,cm,sm,rmm,rmin,rmax,rrr,ar,ac,idbg,nmd=len(cm),knpar=len(rrr))
  
	rrr,ar,ac = sizedis2(kn,cm,sm,rmm,rmin,rmax,nmd=len(cm))


### COMMON blocks:

  /us1/ us11(181,41)
	
  /us2/ us12(181,41)
  
	/us3/ us22(181,41)
  
	/us4/ us33(181,41)
  
	/us5/ us34(181,41)
  
	/us6/ us44(181,41)
  
	/us0/ usea(2,41)
	
### Fortran 90/95 modules:

  mo_par_dls - kn1par,kr1par,km1par,krepar,kimpar,knpar,krpar,kmpar,kmd  
	
	mo_dls - key,key_rd,keyel,keysub,keyls,key_org,key_fx,key_grid1,key_rd1,kn,km,kr,nratn,ndp,wl,rn,rk,pomin,pomax,xext,xabs,xsca,albedo,r,grid,sd,rd,f11,f12,f22,f33,f34,f44,angle,xblr,xldr,distname_o,distname_f,distname_n,comm_name,key_sd,id,nmd,nsd,ksd,cm,sm,rmm,rrr,ar,xgrid,ac  
	
	alloc1 - u11,u12,u22,u33,u34,u44,uea  
	
	alloc - uf11,uf12,uf22,uf33,uf34,uf44,ufea  
	
	mo_intrpl_linear - linear()  
	
	mo_intrpl_spline - intrpl_spline(),e01baf(),e02baf(),e02bbf(),p01abf(),x04aaf(),x04baf(),p01abz()  
	
	phase_func - sint(),sint_spline(),asypar().


# Отрисовка
Отрисовка графиков предполагается с использованием программ plotutils или gnuplot.

