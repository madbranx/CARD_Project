import casadi as casADi

class MaterialProperty:
    CONSTANT = 1
    LINEAR = 2
    POLYNOMIAL = 3

    def __init__(self, log, name, value=None, dependency_type=CONSTANT):
        self.log = log
        self.name = name
        self.value = value
        self.dependency_type = dependency_type

    def get_function(self, temperature=None):
        if self.dependency_type == self.CONSTANT:
            expr = self.value
        elif self.dependency_type == self.LINEAR:
            expr = self.value[0] * temperature + self.value[1]
        elif self.dependency_type == self.POLYNOMIAL:
            expr = sum(coef * (temperature ** i) for i, coef in enumerate(self.value))
        else:
            return None # add Log entry
        return casADi.Function(self.name, [temperature], [expr])


class Component:
    DENSITY = 1
    HEAT_CAPACITY = 2
    THERMAL_CONDUCTIVITY = 3
    COLLISION_AREA = 4
    DIFFUSION_VOLUME = 5
    DYNAMIC_VISCOSITY = 6

    def __init__(self, log, name):
        self.log = log
        self.name = name
        self.collision_area = None
        self.diffusion_volume = None
        self.density = None
        self.heat_capacity = None
        self.thermal_conductivity = None
        self.dynamic_viscosity = None


    def add_property(self, property_type, value=None, dependency_type=MaterialProperty.CONSTANT):
        if property_type == Component.DENSITY:
            name = self.name + '_density'
            self.density = MaterialProperty(self.log, name, value, dependency_type)

        elif property_type == Component.HEAT_CAPACITY:
            name = self.name + '_heat_capacity'
            self.heat_capacity = MaterialProperty(self.log, name, value, dependency_type)

        elif property_type == Component.THERMAL_CONDUCTIVITY:
            name = self.name + '_thermal_conductivity'
            self.thermal_conductivity = MaterialProperty(self.log, name, value, dependency_type)

        elif property_type == Component.COLLISION_AREA:
            name = self.name + '_collision_area'
            self.collision_area = MaterialProperty(self.log, name, value, dependency_type)

        elif property_type == Component.DIFFUSION_VOLUME:
            name = self.name + '_diffusion_volume'
            self.diffusion_volume = MaterialProperty(self.log, name, value, dependency_type)

        elif property_type == Component.DYNAMIC_VISCOSITY:
            name = self.name + '_dynamic_viscosity'
            self.dynamic_viscosity = MaterialProperty(self.log, name, value, dependency_type)

        else:
            self.log.addError('Unknown property type {}'.format(property_type), 3)
            return None

        self.log.addEntry("adding property " + name + " as type " + str(dependency_type), 3)


    def get_density(self, temperature):
        return self.density.get_value(temperature)

    def get_heat_capacity(self, temperature):
        return self.heat_capacity.get_value(temperature)

    def get_thermal_conductivity(self, temperature):
        return self.thermal_conductivity.get_value(temperature)

    def get_collision_area(self):
        return self.collision_area.get_value()

    def get_diffusion_volume(self):
        return self.diffusion_volume.get_value()

    def get_dynamic_viscosity(self):
        return self.dynamic_viscosity.get_value()
