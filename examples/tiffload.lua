-- This is a ARTIMAGEN Lua script that generate the ARTIMAGEN logo SEM image
-- by Petr Cizmar @ NIST

--[[
    As this software was developed as part of work done by the United States
    Government, it is not subject to copyright, and is in the public domain.
    Note that according to GNU.org public domain is compatible with GPL.
]]--

math.randomseed(os.time()) -- Initialization of the random number generator with current time as seed

im=aig_load_image("test_in.tiff");

aig_apply_gaussian_psf(im, 1.5,1,30) -- application of the Gaussian blur

-- definition of the drift/vibration.
freqs = 8 -- the vib. function will be composed of 8 sine functions
ampl=0.2 -- amplitude on each function is 0.2 pixels
freq=100 -- maximum frequency is 100, minimum is 0
vibs={} -- empty table for sine functions
for i = 1, freqs, 1 do
  vibs[i]={math.random()*(ampl), math.random()*(ampl), math.random()*(freq), math.random()*(2*math.pi)} -- random generation of parameters
end

aig_apply_vib(im,10000,50,100,0,vibs) -- application of vibrations
aig_apply_noise(im, "poisson", 200) -- application of Poisson noise
aig_apply_noise(im, "gaussian", 0.1)

--gaig_preview(im) -- gAIG specific! preview. If running with artimagenl, comment this line out.
aig_save_image(im, "test_out.tiff","Logo, by Petr Cizmar") --saving of the image to a file
aig_delete_image(im) -- deletion of the image

