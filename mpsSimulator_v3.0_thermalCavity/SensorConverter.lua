-- rewrite Geo.txt

file = io.open("Sensor.txt","r")

if(file==nil) then
	print(" No file \n")
else

	io.input(file)

		lines = {}
		local pat = "(%S+)%s+(%S+)%s+(%S+)%s*"
		for n1, n2, n3 in string.gfind(io.read("*all"), pat) do
			lines[#lines+1] = string.format("%8.6f %8.6f\n", n1, n3)
		end
	io.close(file)

	file = io.open("Sensor.in","w")
	io.output(file)
		for k,v in pairs(lines) do
			io.write(v)
		end
	io.close(file)

end
