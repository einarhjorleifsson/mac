// srest for fitting SR curves (haked from srmsymc.tpl)
GLOBALS_SECTION
  #include <admodel.h>
  const double pi = 3.141592654;
  
DATA_SECTION
  init_int ybeg;
  init_int yend;
  init_int r; //recruitment age
  init_int A; //plusgroup age
  init_int Ropt; //1=Ricker, 2=B-H, else=SHS (smooth hockey-stick)
  init_int penopt; //0=don't constrain SR params, 1=constrain SR params
  init_matrix srdat(ybeg,yend,1,2);

  int i;
  int n;
  number Rav;
  number Bssbav;
  number Bssbmin;
  number Bssbmax;
  number apinit; //1st stock-recruit parameter (initial value)
  number bpinit; //2nd stock-recruit parameter (initial value)
  number gp; //3rd stock-recruit parameter for SHS model only
  vector R(ybeg,yend); //Recruitment
  vector Bssb(ybeg,yend); //SSB

 LOCAL_CALCS
  n = yend-ybeg-r+1;
  R = column(srdat,1);
  Bssb = column(srdat,2);
  Rav = mean(R);
  Bssbav = mean(Bssb);
  R /= Rav;
  Bssb /= Bssbav;
  Bssbmin = min(Bssb);
  Bssbmax = max(Bssb);
  gp = 0.001;
  if (Ropt == 1) {
    apinit = mfexp(1-Bssb(ybeg));
    bpinit = Bssb(ybeg);
  } else if (Ropt == 2) {
    apinit = Bssb(ybeg)/2.;
    bpinit = (Bssb(ybeg)+1)/2.;
  } else {
    apinit = 1./(1.-gp/2.+sqrt(1.+gp*gp/4.));
    bpinit = 1.;
  }
 END_CALCS

INITIALIZATION_SECTION
  ap apinit;
  bp bpinit;

PARAMETER_SECTION
  init_number ap;
  init_number bp;

  sdreport_number a;
  sdreport_number b;
  sdreport_number sigR;
  sdreport_number scor;
  sdreport_number g;

  vector Rh(ybeg+r,yend);
  vector epsR(ybeg+r,yend);
  number pen;
  sdreport_number neglnL;

  objective_function_value nll;

PROCEDURE_SECTION
  g = gp*Bssbav;
  get_likelihood();
  if (sd_phase()) get_SRpar();

FUNCTION get_likelihood
  dvariable penx;
  pen = 0.;
  penx = posfun(ap,1.e-6,pen);
  penx = posfun(bp,1.e-6,pen);
  if (Ropt == 1)
    Rh = ap*elem_prod(Bssb(ybeg,yend-r),mfexp(-bp*(Bssb(ybeg,yend-r)/Bssb(ybeg)-1.))).shift(Rh.indexmin());
  else if (Ropt == 2) {
    Rh = elem_div(Bssb(ybeg,yend-r),bp+ap*(Bssb(ybeg,yend-r)/Bssb(ybeg)-1.)).shift(Rh.indexmin());
    penx = posfun((bp-ap)*Bssb(ybeg)*Bssbav/ap,1.e-6,pen);
  }
  else {
    Rh = ap*(Bssb(ybeg,yend-r)+sqrt(bp*bp+gp*gp/4.)-sqrt(square(Bssb(ybeg,yend-r)-bp)+gp*gp/4.)).shift(Rh.indexmin());
    penx = posfun(bp-Bssbmin,1.e-6,pen);
    penx = posfun(Bssbmax-bp,1.e-6,pen);
  }
  epsR = log(R(ybeg+r,yend))-log(Rh);
  sigR = norm(epsR)/sqrt(double(n));
  neglnL = 0.5*n*log(2*pi*sigR*sigR) + norm2(epsR)/(2*sigR*sigR);
  if (penopt == 1) neglnL += 1000000000.*pen;
  nll = neglnL;

FUNCTION get_SRpar
  scor = epsR(ybeg+r,yend-1)*epsR(ybeg+r+1,yend).shift(epsR.indexmin())
         /(norm(epsR(ybeg+r,yend-1))*norm(epsR(ybeg+r+1,yend)));
  if (Ropt == 1) {
    a = ap*mfexp(bp)*Rav/Bssbav;
    b = bp/(Bssb(ybeg)*Bssbav);
  } else if (Ropt == 2) {
    a = Bssb(ybeg)*Rav/ap;
    b = (bp-ap)*Bssb(ybeg)*Bssbav/ap;
  } else {
    a = ap*Rav/Bssbav;
    b = bp*Bssbav;
  }

