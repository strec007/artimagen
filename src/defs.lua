-- This is a test aig (LUA) script
--

print("Lua interface");
im=aig_new_image(256,256);
aig_save_image(im,"test.tiff","Yo, man!");
