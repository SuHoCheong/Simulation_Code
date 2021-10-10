
class new_particle:

    def __init__(self,x0,y0,type):
        self.x0 = x0
        self.y0 = y0
        self.newx0 = x0
        self.newy0 = y0
        self.type = type

    def setxcoords(self,x0):
        self.x0 = x0 
 
    def setycoords(self, y0):
        self.y0 = y0

    def setnewxcoords(self,x0):
        self.newx0 = x0 

    def setnewycoords(self, y0):
        self.newy0 = y0

    def getxcoords(self):
        return self.x0

    def getycoords(self):
        return self.y0
 
    def getnewxcoords(self):
        return self.newx0

    def getnewycoords(self):
        return self.newy0

    def gettype(self):
        return self.type

