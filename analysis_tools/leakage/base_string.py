from string import Template
base_string = Template("""\
MCNP6 leakage study radius:${core_radius} cm, height:${core_height} cm. 
c Cell Card
1 1 -6.25 -1 -2 3 imp:n=1  ${comm} cylindrical reactor
2 0 2:-3:(1 -2 3) imp:n=0  ${comm} outside world
 
c Surface Card
1${bc} CZ ${core_radius}
2${bc} PZ ${half_core_height}
3${bc} PZ -${half_core_height}

c Data Card
${fuel_mat}\
kcode ${n_per_cycle} 1 ${non_active_cycles} ${total_cycles}
ksrc 0 0 0
     1 1 1
     -1 -1 -1 
     1 -1 1
     -1 1 1
     1 1 -1
mode n
print
""")

