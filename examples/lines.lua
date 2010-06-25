-- This is a ARTIMAGEN Lua script that generates a set of vertical lines
-- by Petr Cizmar @ NIST

--[[
    As this software was developed as part of work done by the United States
    Government, it is not subject to copyright, and is in the public domain.
    Note that according to GNU.org public domain is compatible with GPL.
]]--

math.randomseed(os.time()) -- Initialization of the random number generator with current time as seed

imsize = 2048 -- Size of a side of the square image
line_width = 20 -- Width of the defined lines
line_pitch = 70 -- Pitch of the defined lines
im = aig_new_image(imsize, imsize) -- Creation of a new empty

background = {} -- New empty table for the background
background_density = 8 -- background "map" will be 8x8 pixels large
for i = 1, background_density*background_density, 1 do
  table.insert(background, math.random()*(0.2-0.1) + 0.1) -- value for each pixel is randomly generated, value varies from 0.1 to 0.2
end
aig_apply_background_image(im, {background_density, background_density}, background) -- application of the background to the image

edge_effect = aig_new_effect("edge", 0.2, 0.3) -- definition of the edge-effect
fine_structure = aig_new_effect("finestructure", 70e-4, 6, 10, 0.95, 1.05) -- definition of the fine structure

features = {} -- new empty table for the feature
i = 0 -- line index definition
while i*line_pitch < imsize do
   curves = {} -- new empty table for crves 
   curves[1] = aig_new_curve("segment", {{-line_width/2,-10}, {-line_width/2,imsize+10}}) -- definition of all curves
   curves[2] = aig_new_curve("segment", {{-line_width/2,imsize+10}, {line_width/2,imsize+10}})
   curves[3] = aig_new_curve("segment", {{line_width/2,imsize+10}, {line_width/2,-10}})
   curves[4] = aig_new_curve("segment", {{line_width/2,-10}, {-line_width/2,-10}})

   lines_feature = aig_new_feature(curves, {edge_effect, fine_structure}, 0.3) -- composition of the curves and effect into a feature
   aig_move_feature(lines_feature, {i*line_pitch,0}) -- shifting of the feature to the center of the image
   table.insert(features, lines_feature); -- insertion of the new feature into the table
   i = i+1
end

lines_sample = aig_new_sample(imsize, imsize, features) -- creation if the sample
aig_paint_sample(im, lines_sample) -- painting of the sample to the image
aig_delete_sample(lines_sample) -- sample is no more needed, thus it should be deleted
aig_apply_gaussian_psf(im, 0.5,1,30) -- application of the Gaussian blur

-- definition of the drift/vibration.
freqs = 8 -- the vib. function will be composed of 8 sine functions
ampl=0.2 -- amplitude on each function is 0.2 pixels
freq=100 -- maximum frequency is 100, minimum is 0
vibs={} -- empty table for sine functions
for i = 1, freqs, 1 do
  vibs[i]={math.random()*(ampl), math.random()*(ampl), math.random()*(freq), math.random()*(2*math.pi)} -- random generation of parameters
end

aig_apply_vib(im,1000,50,100,0,vibs) -- application of vibrations
aig_apply_noise(im, "poisson", 20) -- application of Poisson noise
--aig_apply_noise(im, "gaussian", 0.01)

previm = aig_copy_image(im) -- copy of the original large image
aig_crop_image(previm,0,0,512,512) -- crop it to 512x512 (preview requierement)
gaig_preview(previm) -- gAIG specific - preview. If running with artimagenl, comment this line out.
aig_delete_image(previm) -- deletion of the preview image; it is no more needed

aig_save_image(im, "lines.tiff","Lines, by Petr Cizmar") --saving of the image to a file
aig_delete_image(im) -- deletion of the image