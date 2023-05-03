clear
ifns1 = glob('/lab/tomina/170428/*/*/*coh*.mat');
ifns2 = glob('/lab/tomina/170428/*/*/*/*coh*.mat');
ifns = { ifns1{:}, ifns2{:}}';

F = length(ifns);
for f=1:F
  ifn = ifns{f};
  printf('Loading %i/%i: %s\n', f, F, ifn);
  mc = load(ifn);
  imgs = {coht.extra.img, cohb.extra.img};
  topbot = 1 + (isnan(coht.extra.sig0(300,:)));
  rois = coht.extra.rois;
  sig0 = coht.extra.sig0;
  sig0(:,topbot>1) = cohb.extra.sig0(:,topbot>1);
  ref0 = coht.extra.ref0;
  tt0 = coht.extra.tt0;
  
  expt = mc.x.info.expt;
  tri = mc.x.info.trial;
  printf('Processing %s: %i\n', expt, tri);
  v_ephys = mc.x.analog.dat;
  t_ephys = vscope_ephystime(mc.x);
  ch_ephys = mc.x.analog.info.chanid;
  t_vsd = vscope_ccdtime(mc.x);
  F_vsd1 = mc.roidata{tri};
  F_vsd = F_vsd1(:,:,1);
  isbot = isnan(F_vsd(1,:));
  dFF_vsd(:,isbot) = F_vsd1(:,isbot,2);
  id_vsd = vscope_roiid(1:size(F_vsd,2));
  save(sprintf('/tmp/v2m-%i.mat', tri), '-mat', ...
       't_ephys', 'v_ephys', 'ch_ephys', ...
       't_vsd', 'F_vsd', 'id_vsd', ...
       'expt', 'tri');
end
