
db_ones = {
    "1":  "One",
    "2":  "Two",
    "3":  "Three",
    "4":  "Four",
    "5":  "Five",
    "6":  "Six",
    "7":  "Seven",
    "8":  "Eight",
    "9":  "Nine",
    }
db_teens = {
    "10":  "Ten",
    "11":  "Eleven",
    "12":  "Twelve",
    "13":  "Thirteen",
    "14":  "Fourteen",
    "15":  "Fifteen",
    "16":  "Sixteen",
    "17":  "Seventeen",
    "18":  "Eighteen",
    "19":  "Nineteen",
    }
db_tens = {
    "2":  "Twenty",
    "3":  "Thirty",
    "4":  "Forty",
    "5":  "Fifty",
    "6":  "Sixty",
    "7":  "Seventy",
    "8":  "Eighty",
    "9":  "Ninety",
    }
db_hundreds = {
    "1":  "One Hundred",
    "2":  "Two Hundred",
    "3":  "Three Hundred",
    "4":  "Four Hundred",
    "5":  "Five Hundred",
    "6":  "Six Hundred",
    "7":  "Seven Hundred",
    "8":  "Eight Hundred",
    "9":  "Nine Hundred",
    }
db_magnitudes = {
    "4": "Thousand",
    "7": "Million",
    "10": "Billion",
    "13": "Trillion",
    "16": "Quadrillion",
    "19": "Quintillion",
    "22": "Sextillion",
    "25": "Septillion",
    "28": "Octillion",
    "31": "Nonillion",
    "34": "Decillion",
    "37": "Undecillion",
    "40": "Duodecillion",
    "43": "Tredecillion",
    "46": "Quattuordecillion",
    "49": "Quindecillion",
    "52": "Sexdecillion",
    "55": "Septendecillion",
    "58": "Octodecillion",
    "61": "Novemdecillion",
    "64": "Vigintillion",
    "67": "Unvigintillion",
    "70": "Duovigintillion",
    "73": "Tresvigintillion",
    "76": "Quattuorvigintillion",
    "79": "Quinvigintillion",
    "82": "Sexvigintillion",
    "85": "Septenvigintillion",
    "88": "Octovigintillion",
    "91": "Novemvigintillion",
    "94": "Trigintillion",
    "97": "Untrigintillion",
    "100": "Duotrigintillion",
    "103": "Trestrigintillion",
    "106": "Quattuortrigintillion",
    "109": "Quintrigintillion",
    "112": "Sextrigintillion",
    "115": "Septentrigintillion",
    "118": "Octotrigintillion",
    "121": "Novemtrigintillion",
    "124": "Quadragintillion",
    "127": "Unquadragintillion",
    "130": "Duoquadragintillion",
    "133": "Tresquadragintillion",
    "136": "Quattuorquadragintillion",
    "139": "Quinquadragintillion",
    "142": "Sexquadragintillion",
    "145": "Septenquadragintillion",
    "148": "Octoquadragintillion",
    "151": "Novemquadragintillion",
    "154": "Quinquagintillion",
    "157": "Unquinquagintillion",
    "160": "Duoquinquagintillion",
    "163": "Tresquinquagintillion",
    "166": "Quattuorquinquagintillion",
    "169": "Quinquinquagintillion",
    "172": "Sexquinquagintillion",
    "175": "Septenquinquagintillion",
    "178": "Octoquinquagintillion",
    "181": "Novemquinquagintillion",
    "184": "Sexagintillion",
    "187": "Unsexagintillion",
    "190": "Duosexagintillion",
    "193": "Tresexagintillion",
    "196": "Quattuorsexagintillion",
    "199": "Quinsexagintillion",
    "202": "Sexsexagintillion",
    "205": "Septensexagintillion",
    "208": "Octosexagintillion",
    "211": "Novemsexagintillion",
    "214": "Septuagintillion",
    "217": "Unseptuagintillion",
    "220": "Duoseptuagintillion",
    "223": "Treseptuagintillion",
    "226": "Quattuorseptuagintillion",
    "229": "Quinseptuagintillion",
    "232": "Sexseptuagintillion",
    "235": "Septenseptuagintillion",
    "238": "Octoseptuagintillion",
    "241": "Novemseptuagintillion",
    "244": "Octogintillion",
    "247": "Unoctogintillion",
    "250": "Duooctogintillion",
    "253": "Tresoctogintillion",
    "256": "Quattuoroctogintillion",
    "259": "Quinoctogintillion",
    "262": "Sexoctogintillion",
    "265": "Septenoctogintillion",
    "268": "Octooctogintillion",
    "271": "Novemoctogintillion",
    "274": "Nonagintillion",
    "277": "Unnonagintillion",
    "280": "Duononagintillion",
    "283": "Trenonagintillion",
    "286": "Quattuornonagintillion",
    "289": "Quinnonagintillion",
    "292": "Sexnonagintillion",
    "295": "Septennonagintillion",
    "298": "Octononagintillion",
    "301": "Novemnonagintillion"
}

def main():
    number = input("What's the Number ? ").lower()
    if "e" in number:
        if "." in number:
            float_parser(number)
        else:
            e_parser(number)
    else:
        print(mag_caller(number))

def float_parser(number):
    root, power = number.split("e")
    floater = root.split(".")[1]
    root = root.replace(".", "")
    root = int(root)
    power = int(power)
    power = int(power) - len(floater)
    result = 1
    if power > 0:
        for _ in range(power):
            result = result * 10
        result = result * root
    elif power == 0:
        result = root
    else:
        power = power * -1
        root = str(root)[:-power]
        result = root
    print(f"The Number is : {result:,}")
    print(mag_caller(str(result)))

def e_parser(number):
    root, power = number.split("e")
    root = int(root)
    power = int(power)
    result = 1
    for _ in range(power):
        result = result * 10
    result = result * root
    print(f"The Number is : {result:,}")
    print(mag_caller(str(result)))
    
def mag_caller(number):
    number = number.lstrip("0")
    mag_key = list(db_magnitudes.keys())
    if not number == "":
        if int(number) < 10:
            return ones(number)
        elif 10 <= int(number) < 20:
            return teens(number)
        elif 20 <= int(number) < 100:
            return tens(number)
        elif 100 <= int(number) < 1000:
            return hundreds(number)
        else:
            if len(number) % 3 == 1 or len(number) % 3 == 2:
                mag = (len(number) / 3) - 1
            elif len(number) % 3 == 0:
                mag = (len(number) / 3) - 2

            mag_value = db_magnitudes[mag_key[int(mag)]]
            return computer(number, mag_value, mag_key[int(mag)])
    
def ones(number):
    return db_ones[number]

def teens(number):
    return db_teens[number]

def tens(number):
    for i in db_tens:
        if number[0] == i and number[1] == "0":
            return db_tens[i]
        elif number[0] == i:
            result = db_tens[i]
            return result + " " + str(mag_caller(number[1:]))
            
def hundreds(number):
    for i in db_hundreds:
        if number[0] == i and number[1:len(number)] == "00":
            return db_hundreds[i]
        elif number[0] == i:
            result = db_hundreds[i]
            return result + " " + str(mag_caller(number[1:]))

def computer(number, mag_value, startlen):
    if len(number) == int(startlen):
        for i in db_ones:
            if number[0] == i and set(number[1:]) == {"0"}:
                return f"{db_ones[i]} {mag_value}"
            elif number[0] == i:
                result = f"{db_ones[i]} {mag_value} /"
                return result + " " + str(mag_caller(number[1:]))
    elif len(number) == (int(startlen) + 1):
        if 10 <= int(number[:2]) < 20:
            for i in db_teens:
                if number[:2] == i and set(number[2:]) == {"0"}:
                    return f"{db_teens[i]} {mag_value}"
                elif number[:2] == i:
                    result = f"{db_teens[i]} {mag_value} /"
                    return result + " " + str(mag_caller(number[2:]))
        else:
            for i in db_tens:
                if number[0] == i:
                    if number[1] == "0":
                        if set(number[1:]) == {"0"}:
                            return f"{db_tens[i]} {mag_value}"
                        else:
                            result = f"{db_tens[i]} {mag_value} /"
                            return result + " " + str(mag_caller(number[1:]))
                    else:
                        result = f"{db_tens[i]}"
                        return result + " " + str(mag_caller(number[1:]))
    else:
        for i in db_hundreds:
            if number[0] == i:
                if number[1] == "0":
                    if number[2] == "0":
                        if set(number[2:]) == {"0"}:
                            return f"{db_hundreds[i]} {mag_value}"
                        else:
                            result = f"{db_hundreds[i]} {mag_value} /"
                            return result + " " + str(mag_caller(number[1:]))
                    else:
                        result = f"{db_hundreds[i]}"
                        return result + " " + str(mag_caller(number[1:]))
                elif number[1] != "0":
                    if number[2] == "0":
                        if set(number[2:]) == {"0"}:
                            result = f"{db_hundreds[i]}"
                            return result + " " + str(mag_caller(number[1:]))
                        else:
                            result = f"{db_hundreds[i]}"
                            return result + " " + str(mag_caller(number[1:]))
                    else:
                        result = f"{db_hundreds[i]}"
                        return result + " " + str(mag_caller(number[1:]))

if __name__ == "__main__":
    try:
        main()
    except EOFError:
        pass