-- write Sensor.txt

	dp = 0.01;
	file = io.open("Sensor.in","w")
	io.output(file)
		for y = 0, 1+dp, dp do
			io.write(string.format("%8.6f %8.6f\n", 0.5, y))
		end
		for x = 0, 1+dp, dp do
			io.write(string.format("%8.6f %8.6f\n", x, 0.5))
		end
	io.close(file)
