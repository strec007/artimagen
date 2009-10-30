-- This is a test aig (LUA) script
--

print("Lua interface")
im = aig_new_image(256,256)

bg = {}
dens = 4
for i = 1, dens*dens, 1 do
  table.insert(bg,math.random()*(0.3-0.1) + 0.1)
end

aig_apply_background_image(im,{dens,dens},bg)

ee = aig_new_effect("edge", 0.2, 0.3)

curves = {}
rawset(curves, 1, aig_new_curve("segment",{{0,0},{0,100}}))
rawset(curves, 2, aig_new_curve("segment",{{0,100},{100,100}}))
rawset(curves, 3, aig_new_curve("bezier",{{100,100},{130,80},{130,20},{100,0}}))
rawset(curves, 4, aig_new_curve("segment",{{100,0},{0,0}}))

fs = aig_new_effect("finestructure", 70e-4, 6, 10, 0.95, 1.05)

fe = aig_new_feature(curves, {ee,fs}, 0.3)
aig_move_feature(fe, {100,100})

curves2 = {}
curves2[1] = aig_new_curve("segment", {{0,0}, {50,0}})
curves2[2] = aig_new_curve("segment", {{50,0}, {25,25}})
curves2[3] = aig_new_curve("segment", {{25,25}, {0,0}})
fe2  = aig_new_feature(curves2, {ee, fs}, 0.2)
aig_move_feature(fe2, {50,50})


features = {}
features[1] = fe;
features[2] = fe2;

sa = aig_new_sample(256,256, features)
aig_paint_sample(im, sa)

aig_apply_gaussian_psf(im, 1,1,30)

aig_apply_vib(im,50000,0,0,0,{{1,0.5,1000,0}, {1,2,60,100,50}, {0.3,0.8,5,80,50}})

aig_save_image(im,"test.tiff","Yo, man!")
aig_delete_image(im)
