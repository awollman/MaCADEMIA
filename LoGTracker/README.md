# LoGTracker

An alternative single-channel tracking and stoichiometry pipeline. Detection uses Laplacian
of Gaussian (LoG) filtering by default, and trajectories are built with a Linear Assignment
Problem (LAP) solver with gap closing, rather than the frame-to-frame greedy linking used by
`tracker`. Stoichiometry is calculated by extrapolating trajectory intensity back to t=0
rather than the counting-based approach elsewhere in the repo.

This pipeline does not support dual-channel/ALEX acquisition or cursor mode, and does not
include the rotated-Gaussian (astigmatism-based) PSF fitting available in `tracker` via
`p.GaussSwitch`. Use `tracker` in the Tracking folder if you need those.

## Usage

```
p = createPLoG;
[Spots, frame_average, p, meta_data, image_data] = trackerLoG(fileName, p);

mask = createMask(maskFilePath);
params = createParamsLoG;
[meanStoich, medianStoich, modeStoich, SD, SEM, stoich_all, iSingle] = stoichAnalyser(Spots, mask, params);
```

## trackerLoG

`[Spots, frame_average, p, meta_data, image_data] = trackerLoG(fileName, p)`

Main tracking function. `fileName` is an image stack, folder of tiffs, or an already-loaded
image array. `p` is the parameter structure, created with `createPLoG` if not supplied.
Detection method is set by `p.DetectionMethod`, either `'LoG'` (default) or `'original'`
(uses `findSpots2LoG` and `findSpots3LoG` together). Localisation is set by
`p.LocalisationMethod`, either `'iterative'` (default) or `'quadratic'`. Linking is set by
`p.LinkingMethod`, either `'LAP'` (default, with gap closing via `p.gap_frames`) or `'greedy'`.

Returns `Spots`, an array in the same column format described in the Tracking readme.

## createPLoG

Script that sets the default parameter structure `p` used by `trackerLoG`.

## createParamsLoG

Script that sets the default parameter structure `params` used by `stoichAnalyser` and
`generalLinearStoich`.

## createMask

`mask = createMask(filePath)`

Converts a TIFF image into a binary mask for restricting stoichiometry analysis to a
region of interest.

## generalLinearStoich

`[stoichiometry, iSingle, stoich_all] = generalLinearStoich(Spots, varargin)`

Calculates stoichiometry by fitting a line to the first few frames of each trajectory and
extrapolating intensity back to t=0, then dividing by the single-molecule intensity
`iSingle`. If `iSingle` is not supplied it is estimated from trajectories in the later part
of the movie, where blinking/bleaching means most surviving spots are single molecules.

## stoichAnalyser

`[meanStoich, medianStoich, modeStoich, SD, SEM, stoich_all, iSingle] = stoichAnalyser(Spots, mask, params)`

Wraps `generalLinearStoich`, optionally restricting to trajectories starting inside `mask`,
and reports summary statistics plus a kernel density estimate plot of the stoichiometry
distribution.

## findSpots2LoG / findSpots3LoG

`[result] = findSpots2LoG(frame, method, disk_radius, gaussian, output)`
`[result] = findSpots3LoG(frame, method, disk_radius, gaussian, output)`

Candidate spot detection used by `trackerLoG` when `p.DetectionMethod = 'original'`.
Functionally equivalent to `findSpots2` and `findSpots3` in the Tracking folder, kept as
separate copies here so this pipeline does not depend on files shared with `tracker`.
