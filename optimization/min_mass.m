function test()
mass = test_mass(1100, 1.7e7, 481000, 0.75, 131000);
disp(num2str(mass))
end

function mass = test_mass(T, P, dP, mdot, Q_therm)

mass = py.min_reactor_mass.get_min_mass(T, P, dP, mdot, Q_therm);

end
