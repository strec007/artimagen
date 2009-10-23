-- This is a test aig (LUA) script
--

print("Lua interface")
im = aig_new_image(256,256)

print("Curves")
curves = {}
rawset(curves, 1, aig_new_curve("segment",{{0,0},{0,100}}))
rawset(curves, 2, aig_new_curve("segment",{{0,100},{100,100}}))
rawset(curves, 3, aig_new_curve("bezier",{{100,100},{130,80},{130,20},{100,0}}))
rawset(curves, 4, aig_new_curve("segment",{{100,0},{0,0}}))

print("Feature")
--fe = aig_new_feature({curves[1],curves[2],curves[3],curves[4]});
fe = aig_new_feature(curves);
aig_move_feature(fe,{100,100});
print("Painting feature")
aig_paint_feature(im, fe)
print("Saving image")
aig_save_image(im,"test.tiff","Yo, man!")
aig_delete_image(im)
