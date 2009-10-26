-- This is a test aig (LUA) script
--

print("Lua interface")
im = aig_new_image(256,256)

curves = {}
rawset(curves, 1, aig_new_curve("segment",{{0,0},{0,100}}))
rawset(curves, 2, aig_new_curve("segment",{{0,100},{100,100}}))
rawset(curves, 3, aig_new_curve("bezier",{{100,100},{130,80},{130,20},{100,0}}))
rawset(curves, 4, aig_new_curve("segment",{{100,0},{0,0}}))

fe = aig_new_feature(curves)
aig_move_feature(fe,{100,100})
aig_paint_feature(im, fe)

curves2 = {}
curves2[1] = aig_new_curve("segment",{{0,0},{50,0}})
curves2[2] = aig_new_curve("segment",{{50,0},{25,25}})
curves2[3] = aig_new_curve("segment",{{25,25},{0,0}})
fe2  = aig_new_feature(curves2)
aig_move_feature(fe2, {50,50})
aig_paint_feature(im, fe2)

features = {}
features[1] = fe;
features[2] = fe2;

sa = aig_new_sample(256,256, features)

aig_save_image(im,"test.tiff","Yo, man!")
aig_delete_image(im)
