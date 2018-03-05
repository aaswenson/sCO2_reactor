function mass = calculate_reactor_mass(T, P, dP, m_dot, Q_therm)
    % Call the Python code to calculate the mass of the reactor
    mass = py.min_reactor_mass.get_min_mass(T, P, dP, m_dot, Q_therm);
end