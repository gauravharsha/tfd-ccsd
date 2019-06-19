
class CSStudent: 
    stream = 'cse'     # Class Variable  
    def __init__(self, name, roll): 
        self.name = name  
        self.roll = roll 

    def change_name(self,name):
        self.name = name

a = CSStudent('Gaurav', 1)
print(a.name)
a.change_name('Harsha')
print(a.name)
