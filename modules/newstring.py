def replace(s, old, new):
	return s.replace(old, new)
	
def strip(s, *chars):
	try:
		return s.strip(chars)
	except TypeError:
		return s.strip(None)