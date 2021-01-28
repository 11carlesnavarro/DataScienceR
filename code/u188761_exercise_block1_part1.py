
import math as ma

#### Exercise 1 ####

def get_sphere_volume(radius):
    'Calculates the volume of a sphere with a given radius'
    V = (4/3)*ma.pi*radius**3
    return V

#### Exercise 2 ####

def recursive_factorial(n):
    'Calculates the factorial of a given number recursively'
    if n >= 1:
        return n*recursive_factorial(n-1)
    if n == 0:
        return 1

def factorial(n):
    'Calculates the factorial of a given number'
    factorial = 1
    while n > 0:
        factorial = factorial*n
        n = n - 1     
    return factorial
        
#### Exercise 3 ####

def recursive_count_up(n, odd):
    'Count up numbers from 0 to n'
    if n >= 0:
        count = recursive_count_up(n-1, odd)+1
        if odd == True:
            if (count % 2) != 0:
                print (count)
        if odd == False:
            print (count)
    return n

def count_up(n, odd):
    m = 0
    if odd == True:
        while m <= n:
            if (m % 2) != 0:
                print(m)
            m = m+1
    if odd == False:
        while m <= n:
            print (m)
            m = m+1

#### Exercise 4 ####

def get_final_price(price, discount_percentage = 20):                       # I put price on the first place to avoid errors when calling the function.
    """Return the final price after applying the discount percentage"""
    return (price - (price * discount_percentage/100))                      # I change "10" by 10 because an integer it is needed. Instead we could use the function  
                                                                            # int() inside to transform the string into an integer. 
                                                                            # Finally, I change the function calculation in order to return the true final price.

