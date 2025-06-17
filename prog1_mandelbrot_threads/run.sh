echo "view = 1"
for i in {2..8}; do echo "Testing with $i threads:"; ./mandelbrot --threads $i --view 1; done
