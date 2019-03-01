while True:
    try:
        search_field = int(input("Please enter an integer: "))
        search_field = int(search_field)
        break
    except ValueError:
        print("No valid integer! Please try again ...")
print("Great, you successfully entered an integer!")
###########
search_field_no = 0

search_term = input(
    'Property to search for: \nCapitalise for function, for example:\'Biosynthesis\' instead of \'biosynthesis\': ')

file_name = input("Input name of genbank file (make sure it's in the working directory): ")

while search_field_no == 0:
    try:
        search_field_no = int(input('Search in: 1 = function, 2 = product: '))
    except ValueError:
        print('Please choose between option 1 or 2')
    if search_field_no == 1:
        search_field = 'function'
        break
    elif search_field_no == 2:
        search_field = 'product'
        break
    else:
        print('Please type in 1 or 2 to choose.')
        break
return (file_name, search_term, search_field)