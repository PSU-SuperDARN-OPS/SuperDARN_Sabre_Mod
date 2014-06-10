# runs nec2/mp over multiple files in parallel
# ram usage for running 4 at once is ~4 gb, loads ~60% of a i7-4700hq on average
# view out files with xnecview -z0 100 -log {filename}
python2 log_antenna.py
find *.nec | parallel -j 4 'wine nec2mp/nec2dxs11k.exe {} out/{}.out'
