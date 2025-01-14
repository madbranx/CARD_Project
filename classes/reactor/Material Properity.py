class MaterialProperty:
    def __init__(self, value, dependency_type='constant', coefficients=None):
        self.value = value
        self.dependency_type = dependency_type
        self.coefficients = coefficients if coefficients else []

    def get_value(self, temperature):
        if self.dependency_type == 'constant':
            return self.value
        elif self.dependency_type == 'linear':
            return self.coefficients[0] * temperature + self.coefficients[1]
        elif self.dependency_type == 'polynomial':
            return sum(coef * (temperature ** i) for i, coef in enumerate(self.coefficients))
        else:
            raise ValueError("Invalid dependency type")

class Material:
    def __init__(self, name, density, heat_capacity, thermal_conductivity):
        self.name = name
        self.density = MaterialProperty(density)
        self.heat_capacity = MaterialProperty(heat_capacity)
        self.thermal_conductivity = MaterialProperty(thermal_conductivity)

    def add_property(self, property_name, value, dependency_type='constant', coefficients=None):
        if property_name == 'density':
            self.density = MaterialProperty(value, dependency_type, coefficients)
        elif property_name == 'heat_capacity':
            self.heat_capacity = MaterialProperty(value, dependency_type, coefficients)
        elif property_name == 'thermal_conductivity':
            self.thermal_conductivity = MaterialProperty(value, dependency_type, coefficients)
        else:
            raise ValueError("Invalid property name")

    def get_property(self, property_name, temperature):
        if property_name == 'density':
            return self.density.get_value(temperature)
        elif property_name == 'heat_capacity':
            return self.heat_capacity.get_value(temperature)
        elif property_name == 'thermal_conductivity':
            return self.thermal_conductivity.get_value(temperature)
        else:
            raise ValueError("Invalid property name")

# Example usage:
material = Material(7800, 500, 45)
material.add_property('density', 7800, 'linear', [0.1, 7800])
material.add_property('heat_capacity', 500, 'polynomial', [0.01, 0.1, 500])