from location_objects import (SingleBaseLocation, RangeLocation,
							  RemoteLocation, AdjoiningLocation,
							  JoinedLocation, ComplementLocation)


def parse_location(location_string):
	location = __parse_string_arguments(location_string.strip())
	if len(location) > 1:
		raise ValueError('Cannot parse {} to a Location!'
						 .format(location_string))
	# Check for a function
	result = __execute_function(location[0])
	return result or __parse_location_string(location[0])


def __execute_function(argument, top_level=True):
	for func_map in function_mapping:
		if argument.startswith(func_map['name']):
			remaining = argument[len(func_map['name']):].lstrip()
			if remaining[0] != '(':
				raise ValueError('Expected \'()\'')
			remaining = remaining[1:]
			arguments = __parse_string_arguments(remaining)
			if func_map['one_arg'] and len(arguments) > 1:
				raise ValueError('1 argument allowed for ' + func_map['name'])
			if top_level:
				for i in range(len(arguments)):
					result = __execute_function(arguments[i], False)
					if result is not None:
						arguments[i] = result
			for i in range(len(arguments)):
				if isinstance(arguments[i], str):
					arguments[i] = __parse_location_string(arguments[i])
			return func_map['function'](*arguments)


def __parse_string_arguments(location_string):
	argument = ""
	arguments = []
	parentheses_level = 0
	for char in location_string:
		if char == '(':
			parentheses_level += 1
			argument += char
		elif char == ')' and parentheses_level != 0:
			parentheses_level -= 1
			argument += char
		elif char == ')' and parentheses_level == 0:
			break
		elif char == ',' and parentheses_level == 0:
			arguments.append(argument.strip())
			argument = ""
		else:
			argument += char
	arguments.append(argument.strip())
	return arguments


def __parse_location_string(location_string):
	for (char, cls) in location_mapping:
		if char in location_string:
			return cls(location_string)
	return SingleBaseLocation(location_string)


location_mapping = ((':', RemoteLocation), ('..', RangeLocation),
					('^', AdjoiningLocation))


def __create_func_mapping(name, func, one_arg=False, sub_func=False,
						  contain_sub_func=False):
	return dict(name=name, function=func, one_arg=one_arg, sub_func=sub_func,
				contain_sub_func=contain_sub_func)


function_mapping = [
	__create_func_mapping('complement', ComplementLocation, True, False,
						  True),
	__create_func_mapping('join', JoinedLocation, False, True, False),
	# __create_func_mapping('order', __parse_order, False, True, False)
]
