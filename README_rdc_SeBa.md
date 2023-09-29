**Reducing SeBa data: rdc_SeBa**

The program rdc_SeBa reduces the output file from SeBa. You can find it in SeBa's directory rdc. One of the big advantages of using rdc_SeBa is that it picks out only those binary systems that you're interested in, and reduces the output to a regular number of lines, which makes the data handeling easier. 

To run: 
```
less SeBa.data | ./rdc_SeBa -f -p white_dwarf -s white_dwarf 
```
or
```
less SeBa.data | ./rdc_SeBa -f -p white_dwarf -s white_dwarf > SeBa_wdwd.data
```

*Optional parameters:*
```
-f  first occasion - for a given binary system, the first occasion that fits the parameters. In case of double white dwarfs, this will be the formation of the double white dwarf system, and later stages of the double white dwarf system will not be shown. 
-R  full evolution - for a given binary that fulfills the criteria, every line in SeBa.data is printed.     
-p  stellar type of primary star
-s  stellar type of secondary star
-B binary_type
```

Options for the stellar types are: 
```
any, proto_star, planet, brown_dwarf, main_sequence, hertzsprung_gap, sub_giant, horizontal_branch, super_giant, helium_star, helium_giant, white_dwarf, neutron_star, black_hole
ps, pl, bd, ms, hg, gs, hb, sg, he, gh, hd, cd, od (last three are types of white dwarfs), ns, bh
```

Options for the binary types are: 
```
detached [default], semi_detached, contact, common_envelope, double_spiral_in, merged, disrupted, spiral_in
```
with: semi_detached meaning stable mass transfer, common_envelope meaning gamma common envelope, and spiral_in meaning alpha common envelope.


Other options constraining binary parameters (not used for long time): 
```
-a,-A   min/max separation
-m,-M   min/max primary mass
-n,-N   min/max secondary mass
-q,-Q   min/max mass ratio
-e,-E   min/max eccentricity
```