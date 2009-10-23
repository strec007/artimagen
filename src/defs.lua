-- This is a test aig (LUA) script
--

print("Lua interface")
im = aig_new_image(256,256)

curve1 = aig_new_curve("segment",{{0,0},{0,100}})
curve2 = aig_new_curve("segment",{{0,100},{100,0}})
curve3 = aig_new_curve("segment",{{100,0},{0,0}})

aig_save_image(im,"test.tiff","Yo, man!")
aig_new_feature({curve1,curve2,curve3});
aig_delete_image(im)
