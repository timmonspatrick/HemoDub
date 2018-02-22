import requests
for n in range(1, 10893):
    try:
        url = "https://dbaasp.org/api/v1?query=peptide_card&peptide_id=%s&format=json" % str(n)
        
        data = requests.get(url).text
        with open("DBAASP/DBAASP%s.html" % (str(n)), "w", encoding="utf-8") as f:
            f.write(data)
        print("DBAASP " + str(n) + " successful")
    except:
        print("DBAASP " + str(n) + " failed")

for n in range(1001, 4202):
    try:
        url = "http://crdd.osdd.net/raghava/hemolytik/display.php?details=%s" %(str(n))
        data = requests.get(url)
        data.encoding = "utf-8"
        data = data.text
        with open("HEMOLYTIK/HEMOLYTIK%s.html" % (str(n)), "w", encoding="utf-8") as f:
            f.write(data)
        print("HEMOLYTIK " + str(n) + " successful")
    except:
        print("HEMOLYTIK " + str(n) + " failed")
