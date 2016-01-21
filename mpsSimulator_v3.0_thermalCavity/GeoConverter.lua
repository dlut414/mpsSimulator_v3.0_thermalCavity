-- rewrite Geo.txt

file = io.open("Geo.txt","r")

if(file==nil) then
	print(" No file \n")
else

	io.input(file)
		local L1, L2, L3
		L1 = io.read("*line")
		L2 = io.read("*line")
		L3 = io.read("*line")

		lines = {}
		local pat = "(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s*"
		for n1, n2, n3, n4, n5, n6, n7 in string.gfind(io.read("*all"), pat) do
			lines[#lines+1] = string.format("%2d %8.6f %8.6f %8.6f %8.6f\n", n1, n2, n4, n5, n7)
		end
	io.close(file)

	file = io.open("Geo.in","w")
	io.output(file)
		io.write(L1,"\n")
		io.write(L2,"\n")
		io.write(L3,"\n")
		for k,v in pairs(lines) do
			io.write(v)
		end
	io.close(file)

end
