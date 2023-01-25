import requests
from bs4 import BeautifulSoup


class SpeechModelScrape:
    def __init__(self):
        self.url = "https://alphacephei.com/vosk/models"
        self.page_html = requests.get(self.url).text

        self.soup = BeautifulSoup(self.page_html, "lxml")
        # print(self.soup.prettify())  # print the parsed data of html
        # print(self.soup.title)

        self.model_table = self.soup.find("table", attrs={"class": "table table-bordered"})

        self.headings = []
        for th in self.model_table.thead.find_all("th"):
            self.headings.append(th.text.replace("\n", "").strip())
        print(self.headings)

        self.data = []
        l_models = []
        init_flag = 1
        for tr in self.model_table.tbody.find_all("tr"):
            tr_model = []
            c = 0
            for td in tr:  # for column in row
                if c % 2 != 0:  # even rows are dividers and empty
                    tr_model.append(td.text.replace('\n', ' ').strip())
                c += 1
            if tr_model[4] == "Â":  # start of a new language
                if init_flag:
                    l_models.append([tr_model[0]])
                else:
                    self.data.append(l_models.copy())  # current contents of l_models make up one language of models
                    l_models.clear()  # clear l_models for next language
                    l_models.append([tr_model[0]])  # remove blank columns in model language row
            else:
                l_models.append(tr_model)
            init_flag = 0

        # print(self.data)
        # for language in self.data:
            # print(language[0][0])
            # print(language[1:])

        print(self.data[0][1:][0][0])


if __name__ == "__main__":
    SpeechModelScrape()
