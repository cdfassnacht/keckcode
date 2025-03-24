

class LRIS_MOSmask:

    def __init__(self,maskname,RA=None,Dec=None):
        self.observations = {}
        self.maskname = maskname
        self.RA = RA
        self.Dec = Dec
        self.final = None

    def add_observation(self,obsname,observation):
        self.observations[obsname] = observation

    def setRA(self,RA):
        self.RA = RA

    def setDec(self,Dec):
        self.Dec = Dec

    def process_observation(self,obsname):
        self.observations[obsname] = pipeline(self.observations[obsname])

    def process_observations(self):
        self.final = uber_pipeline(self.observations)

    def combine_observations(self):
        self.final = combine_obs(self.observations)

    def publish(self):
        if self.final is None:
            print 'No final products to publish.'
            return

