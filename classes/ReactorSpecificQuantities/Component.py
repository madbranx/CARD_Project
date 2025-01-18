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

class Component:
    DENSITY = 1
    HEAT_CAPACITY = 2
    THERMAL_CONDUCTIVITY = 3
    COLLISION_AREA = 4
    DIFFUSION_VOLUME = 5

    def __init__(self, name):
        self.name = name
        self.collision_area = None
        self.diffusion_volume = None
        self.density = None
        self.heat_capacity = None
        self.thermal_conductivity = None

    def add_property(self, property_name, value, dependency_type='constant', coefficients=None):
        if property_name == Component.DENSITY:
            self.density = MaterialProperty(value, dependency_type, coefficients)

        elif property_name == Component.HEAT_CAPACITY:
            self.heat_capacity = MaterialProperty(value, dependency_type, coefficients)

        elif property_name == Component.THERMAL_CONDUCTIVITY:
            self.thermal_conductivity = MaterialProperty(value, dependency_type, coefficients)

        elif property_name == Component.COLLISION_AREA:
            self.collision_area = value

        elif property_name == Component.DIFFUSION_VOLUME:
            self.diffusion_volume = value


    def get_density(self, temperature):
        return self.density.get_value(temperature)

    def get_heat_capacity(self, temperature):
        return self.heat_capacity.get_value(temperature)

    def get_thermal_conductivity(self, temperature):
        return self.thermal_conductivity.get_value(temperature)

    def get_collision_area(self):
        return self.collision_area

    def get_diffusion_volume(self):
        return self.diffusion_volume
