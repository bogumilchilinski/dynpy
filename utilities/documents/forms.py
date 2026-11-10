# from .guides import (UsageOfDynamicSystemsGuide, Guide, EngeneeringDrawingGuide, DevelopmentGuide, IntroToPandasGuide,
#                     BasicsOfODESystemGuide, BasicsOfDynSysImplementationGuide, BasicsOfReportingGuide, ResearchProjectGuidelines,
#                     InterimProjectGuidelines,IntroDynPyProjectGuidelines, BasicsOfReportComponentImplementationGuide,
#                     GithubSynchroGuide, IntroToCocalcGuide)
# from sympy import *
import datetime
import os
import shutil
import csv
from typing import List, Optional, Union
import pandas as pd

from pylatex import (  # Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
    Command,
    Document,
    NewPage,
    Package,
    Tabularx,
)
from pylatex.base_classes import Environment

# from pylatex.section import Paragraph, Chapter
from pylatex.utils import NoEscape  # italic,

from ..report import (
    CurrentContainer,
    IPMarkdown,
    Markdown,
    ObjectCode,
    ReportText,
    display,
)

class ThesisForm(Document):
    import datetime
import os
import shutil
import csv
from typing import List, Optional, Union
import pandas as pd

from pylatex import (  # Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
    Command,
    Document,
    NewPage,
    Package,
    Tabularx,
    Tabular
)
from pylatex.base_classes import Environment

# from pylatex.section import Paragraph, Chapter
from pylatex.utils import NoEscape  # italic,

from dynpy.utilities.report import *

class ThesisForm(Document):

    packages = [
                Package("polski"),
                Package("setspace"),
                Package("array"),
                Package("ragged2e"),
                Package("tgtermes"),
                Package("tabularx"),
                Package("ulem", options="normalem"),
                Package("xcolor"),
                Package("textpos", options=["absolute", "overlay"]),
            ]
        
    
    """Generate a LaTeX thesis form card from a row of thesis 
    ."""

    def __init__(
        self,
        default_filepath: str = "default_filepath",
        *,
        dane: Optional[Union[dict, pd.Series]] = None,
        documentclass: str = "article",
        document_options: Optional[List[str]] = None,
        fontenc: str = "T1",
        inputenc: str = "utf8",
        font_size: str = "normalsize",
        lmodern: bool = False,
        textcomp: bool = True,
        microtype: bool = False,
        page_numbers: bool = False,
        indent: Optional[Union[str, int]] = None,
        geometry_options: Optional[List[str]] = None,
        data: Optional[dict] = None,
        auto_build: bool = True,
        auto_generate: bool = False,
        clean_tex: bool = False,
    ):
        if geometry_options is None:
            geometry_options = [
                "left=2.5cm",
                "right=2.5cm",
                "top=2.5cm",
                "bottom=2.5cm",
            ]

        super().__init__(
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )

        self.base_name: Optional[str] = None
        self.tex_path: Optional[str] = None
        self.pdf_path: Optional[str] = None

        if dane is not None and auto_build:
            self.build(dane)
            if auto_generate:
                self.generate_tex()

    @staticmethod
    def bezpieczny_tex(value: object) -> str:
        if value is None:
            return ""

        try:
            if pd.isna(value):
                return ""
        except Exception:
            pass

        text = str(value).strip()
        if not text:
            return ""

        replacements = {
            "\\": r"\\textbackslash{}",
            "&": r"\\&",
            "%": r"\\%",
            "$": r"\\$",
            "#": r"\\#",
            "_": r"\\_",
            "{": r"\\{",
            "}": r"\\}",
            "~": r"\\textasciitilde{}",
            "^": r"\\^{}",
        }
        for key, replacement in replacements.items():
            text = text.replace(key, replacement)

        return text

    @staticmethod
    def _get_value(dane: object, key: str, default: str = "") -> str:
        if dane is None:
            return default

        if isinstance(dane, dict):
            return dane.get(key, default)

        try:
            if key in dane:
                return dane[key]
        except Exception:
            pass

        return default

    @staticmethod
    def _get_indexed_value(dane: object, index: int, default: str = "") -> str:
        try:
            if hasattr(dane, "iloc"):
                return dane.iloc[index]
        except Exception:
            pass

        return default

    def build(self, dane: Union[dict, pd.Series]) -> Optional[str]:
        dyplomant = self.bezpieczny_tex(self._get_value(dane, "Dyplomant"))
        album = self.bezpieczny_tex(self._get_value(dane, "Numer albumu"))

        if not dyplomant or not album:
            dyplomant = 'BCh'
            album ='1234'

        print(f"\n--- Generowanie bazy dla: {dyplomant} ({album}) ---")
        self.base_name = f"karta_{album}_{dyplomant.replace(' ', '_')}"
        self.tex_path = f"./output/{self.base_name}"
        self.pdf_path = f"{self.tex_path}.pdf"

        kierunek = self.bezpieczny_tex(self._get_value(dane, "Kierunek"))
        specjalnosc = self.bezpieczny_tex(self._get_value(dane, "Specjalność"))
        email = self.bezpieczny_tex(self._get_value(dane, "Adres e-mail"))
        prowadzacy = self.bezpieczny_tex(self._get_value(dane, "Prowadzący"))
        konsultant = self.bezpieczny_tex(self._get_value(dane, "Konsultant"))
        temat_tresc = self.bezpieczny_tex(self._get_value(dane, "Treść tematu"))
        opis = self.bezpieczny_tex(self._get_value(dane, "Opis")).replace("\n", " ")
        jezyk = self.bezpieczny_tex(self._get_value(dane, "Język"))
        data_wydania = self.bezpieczny_tex(self._get_value(dane, "Data wydania"))

        rodzaj_studiow = str(self._get_indexed_value(dane, 1)).lower()
        if "niestacjonarne" in rodzaj_studiow:
            tex_rodzaj = r"\textbf{\sout{Stacjonarne} / Niestacjonarne (zaoczne) *}"
        else:
            tex_rodzaj = r"\textbf{Stacjonarne / \sout{Niestacjonarne (zaoczne)} *}"

        stopien_studiow = str(self._get_indexed_value(dane, 2)).lower()
        if "magisterskie" in stopien_studiow or "ii" in stopien_studiow:
            tex_stopien = r"\textbf{\sout{I (inżynierskie)} / II (magisterskie) *}"
        else:
            tex_stopien = r"\textbf{I (inżynierskie) / \sout{II (magisterskie)} *}"

        status_tematu = str(self._get_indexed_value(dane, 10)).lower()
        if "zmieniony" in status_tematu:
            tex_status = r"\noindent Temat pracy dyplomowej (\sout{nowy} / zmieniony)*:\\[0.1cm]"
        else:
            tex_status = r"\noindent Temat pracy dyplomowej (nowy / \sout{zmieniony})*:\\[0.1cm]"

        zadania_kody = []
        for i in range(1, 6):
            zad_kol = f"Zadanie {i}"
            if zad_kol in dane and pd.notna(dane[zad_kol]) and str(dane[zad_kol]).strip():
                zad_tekst = self.bezpieczny_tex(dane[zad_kol])
                zadania_kody.append(f"{i}. {zad_tekst}" + r" \dotfill\\[0.1cm]")

        if zadania_kody:
            zadania_kody[-1] = zadania_kody[-1].replace(r"\\[0.1cm]", r"\\[0.3cm]")

        self.append(NoEscape(r"\setstretch{1.15}"))
        self.append(NoEscape(r"\pagestyle{empty}"))

        self.append(NoEscape(r"\noindent\textcolor{gray}{Wydział Samochodów i Maszyn Roboczych}\par\vspace{0.5cm}"))
        self.append(NoEscape(r"\begin{center} \Large \textbf{Karta pracy dyplomowej} \end{center}"))
        self.append(Command("vspace", "0.5cm"))

        with self.create(Tabular("p{3.5cm} p{11.5cm}")) as table:
            table.add_row(("Rodzaj studiów:", NoEscape(tex_rodzaj)))
            table.add_row(("Stopień studiów:", NoEscape(tex_stopien)))
            table.add_row(("Kierunek:", NoEscape(kierunek + r" \dotfill")))
            table.add_row(("Specjalność:", NoEscape(specjalnosc + r" \dotfill")))
            table.add_row(("Dyplomant:", NoEscape(dyplomant + r" \dotfill")))
            table.add_row(("Numer albumu:", NoEscape(album + r" \dotfill")))
            table.add_row(("Adres e-mail:", NoEscape(email + r" \dotfill")))
            table.add_row(("Prowadzący:", NoEscape(prowadzacy + r" \dotfill" if prowadzacy else r"\dotfill")))
            table.add_row(("Konsultant:", NoEscape(konsultant + r" \dotfill" if konsultant else r"\dotfill")))

        self.append(NoEscape(r"\par\vspace{0.5cm}"))
        self.append(NoEscape(tex_status))
        
        # Zastosowanie bezpiecznego łączenia stringów w treści
        self.append(NoEscape(temat_tresc + r" \dotfill\\[0.3cm]"))
        self.append(NoEscape(r"\noindent Zwięzły opis celu pracy:\\[0.1cm]"))
        self.append(NoEscape(opis + r" \dotfill\\[0.3cm]"))
        
        self.append(NoEscape(r"\noindent Główne zadania do wykonania:\\[0.1cm]"))
        for zad in zadania_kody:
            self.append(NoEscape(zad + "\n"))

        self.append(NoEscape(rf"\noindent Język opracowania: {jezyk}" + r" \dotfill\\[1cm]"))

        self.append(NoEscape(rf'''
        \begin{{textblock*}}{{16cm}}(2.5cm, 22.0cm)
            \begin{{tabularx}}{{16cm}}{{>{{\centering\arraybackslash}}X >{{\centering\arraybackslash}}X}}
                ........................................ & .......... {data_wydania} .......... \\
                {{\footnotesize Podpis studenta}} & {{\footnotesize Data wydania}} \\[1.2cm]
                ........................................ & ........................................ \\
                {{\footnotesize Podpis prowadzącego pracę}} & {{\footnotesize Podpis prodziekana}} \\
            \end{{tabularx}}
        \end{{textblock*}}
        ''') )

        self.append(NoEscape(r"\vfill"))
        self.append(NoEscape(r"{\small *) niepotrzebne skreślić}"))

        self.append(Command("newpage"))
        self.append(NoEscape(r"\noindent\textcolor{gray}{Wydział Samochodów i Maszyn Roboczych}\par\vspace{0.5cm}"))
        self.append(NoEscape(r"\noindent \textbf{Konsultacje i ocena zaawansowania pracy}\\[0.5cm]"))
        self.append(NoEscape(r"\renewcommand{\arraystretch}{2}"))

        with self.create(Tabularx(NoEscape(r'|p{2cm}|X|p{3.5cm}|p{4.5cm}|'))) as consult_table:
            consult_table.add_hline()
            consult_table.add_row(("Data", "Temat konsultacji", "Zaawansowanie pracy", "Podpis prowadzącego lub konsultanta"))
            consult_table.add_hline()
            for _ in range(13):
                consult_table.add_row(("", "", "", ""))
                consult_table.add_hline()

        return self
    
    def sign_it(self,signiture):
        import fitz
        import pandas as pd
        import os

        self.generate_pdf()

        nazwa_bazowa = self.base_name
        pdf_bazowy = self.pdf_path
        pdf_wyjsciowy = f"{self.pdf_path}".replace('.pdf','_bch_signed.pdf')
        podpis_promotor = signiture

        if not os.path.exists(pdf_bazowy):
            print(f"Brak pliku docelowego: {pdf_bazowy}")
            return

        print(f"Dodawanie podpisu promotora do: {self.base_name}")

        try:
            doc_pdf = fitz.open(pdf_bazowy)
            strona = doc_pdf[0]

            prostokat_promotor = fitz.Rect(122, 657, 272, 697)

            if os.path.exists(podpis_promotor):
                strona.insert_image(prostokat_promotor, filename=podpis_promotor)
                doc_pdf.save(pdf_wyjsciowy)
                print(f"  -> Zapisano jako: {pdf_wyjsciowy}")
            else:
                print(f"  -> UWAGA: Brak obrazu {podpis_promotor}! Pominięto.")
                doc_pdf.save(pdf_wyjsciowy) 

            doc_pdf.close()
        except Exception as e:
            print(f"Błąd przy wstawianiu podpisu dla {self.base_name}: {e}")
    
    
import pandas as pd

#Tutaj wpisz swoje dane
dane_dict = {
    "Dyplomant": "Jan Kowalski",
    "Numer albumu": "123456",
    "Kierunek": "Mechanika i Budowa Maszyn",
    "Specjalność": "Pojazdy samochodowe",
    "Adres e-mail": "jan.kowalski@student.edu.pl",
    "Prowadzący": "dr inż. Nowak",
    "Konsultant": "mgr inż. Wiśniewski",
    "Treść tematu": "Analiza wytrzymałościowa ramy pojazdu",
    "Opis": "Celem pracy jest analiza wytrzymałościowa konstrukcji ramy pojazdu z wykorzystaniem MES.",
    "Język": "polski",
    "Data wydania": "01.10.2025",

    "Zadanie 1": "Przegląd literatury",
    "Zadanie 2": "Budowa modelu numerycznego",
    "Zadanie 3": "Przeprowadzenie symulacji",
    "Zadanie 4": "Analiza wyników",
    "Zadanie 5": "Opracowanie wniosków",
}

# część indeksowa (symuluje pd.Series z kolumnami bez nazw)
# indeks 1 -> rodzaj studiów
# indeks 2 -> stopień
# indeks 10 -> status tematu
#dane_lista = [
#    None,
#    "stacjonarne",      # indeks 1
#    "I stopnia",        # indeks 2
#    None, None, None, None, None, None,
#    None,
#    "nowy"              # indeks 10
#]



dane = pd.Series(dane_dict)
# for i, val in enumerate(dane_lista):
#     if val is not None:
#         dane.loc[i] = val


#Wywołanie
#doc = ThesisForm(
#   default_filepath="./output/karta_123456_Jan_Kowalski",
#   dane=dane_dict,
#   auto_build=True,
#   auto_generate=True  # wygeneruje .tex
#)

doc.generate_pdf(clean_tex=False)
    pass
    

class InterimProjectForm(ThesisForm):
    import datetime
import os
import shutil
import csv
from typing import List, Optional, Union
import pandas as pd

from pylatex.utils import bold

from pylatex import (  # Section, Subsection, Subsubsection, Itemize,  HorizontalSpace, Description, Marker
    Command,
    Document,
    NewPage,
    Package,
    Tabularx,
    Tabular
)
from pylatex.base_classes import Environment

# from pylatex.section import Paragraph, Chapter
from pylatex.utils import NoEscape  # italic,

from dynpy.utilities.report import *

class InterimProjectForm(Document):

    packages = [
                Package("polski"),
                Package("setspace"),
                Package("array"),
                Package("ragged2e"),
                Package("tgtermes"),
                Package("tabularx"),
                Package("ulem", options="normalem"),
                Package("xcolor"),
                Package("parskip"),
                Package("textpos", options=["absolute", "overlay"]),
            ]
        
    
    """Generate a LaTeX thesis form card from a row of thesis 
    ."""

    def __init__(
        self,
        default_filepath: str = "default_filepath",
        *,
        dane: Optional[Union[dict, pd.Series]] = None,
        documentclass: str = "article",
        document_options: Optional[List[str]] = None,
        fontenc: str = "T1",
        inputenc: str = "utf8",
        font_size: str = "normalsize",
        lmodern: bool = False,
        textcomp: bool = True,
        microtype: bool = False,
        page_numbers: bool = False,
        indent: Optional[Union[str, int]] = None,
        geometry_options: Optional[List[str]] = None,
        data: Optional[dict] = None,
        auto_build: bool = True,
        auto_generate: bool = False,
        clean_tex: bool = False,
    ):
        if geometry_options is None:
            geometry_options = [
                "left=2.5cm",
                "right=2.5cm",
                "top=2.5cm",
                "bottom=2.5cm",
            ]

        super().__init__(
            default_filepath=default_filepath,
            documentclass=documentclass,
            document_options=document_options,
            fontenc=fontenc,
            inputenc=inputenc,
            font_size=font_size,
            lmodern=lmodern,
            textcomp=textcomp,
            microtype=microtype,
            page_numbers=page_numbers,
            indent=indent,
            geometry_options=geometry_options,
            data=data,
        )

        self.base_name: Optional[str] = None
        self.tex_path: Optional[str] = None
        self.pdf_path: Optional[str] = None

        if dane is not None and auto_build:
            self.build(dane)
            if auto_generate:
                self.generate_tex()

    @staticmethod
    def bezpieczny_tex(value: object) -> str:
        if value is None:
            return ""

        try:
            if pd.isna(value):
                return ""
        except Exception:
            pass

        text = str(value).strip()
        if not text:
            return ""

        replacements = {
            "\\": r"\\textbackslash{}",
            "&": r"\\&",
            "%": r"\\%",
            "$": r"\\$",
            "#": r"\\#",
            "_": r"\\_",
            "{": r"\\{",
            "}": r"\\}",
            "~": r"\\textasciitilde{}",
            "^": r"\\^{}",
        }
        for key, replacement in replacements.items():
            text = text.replace(key, replacement)

        return text

    @staticmethod
    def _get_value(dane: object, key: str, default: str = "") -> str:

        return dane.get(key, default)


    @staticmethod
    def _get_indexed_value(dane: object, index: int, default: str = "") -> str:
        try:
            if hasattr(dane, "iloc"):
                return dane.iloc[index]
        except Exception:
            pass

        return default

    def build(self, dane: Union[dict, pd.Series]) -> Optional[str]:
        dyplomant = self.bezpieczny_tex(self._get_value(dane, "Dyplomant"))
        album = self.bezpieczny_tex(self._get_value(dane, "Numer albumu"))

        if not dyplomant or not album:
            dyplomant = 'BCh'
            album ='1234'

        print(f"\n--- Generowanie bazy dla: {dyplomant} ({album}) ---")
        self.base_name = f"karta_{album}_{dyplomant.replace(' ', '_')}"
        self.tex_path = f"./output/{self.base_name}"
        self.pdf_path = f"{self.tex_path}.pdf"

        self.default_filepath = self.tex_path
        
        kierunek = self.bezpieczny_tex(self._get_value(dane, "Kierunek"))
        specjalnosc = self.bezpieczny_tex(self._get_value(dane, "Specjalność"))
        email = self.bezpieczny_tex(self._get_value(dane, "Adres e-mail"))
        prowadzacy = self.bezpieczny_tex(self._get_value(dane, "Prowadzący"))
        konsultant = self.bezpieczny_tex(self._get_value(dane, "Konsultant"))
        temat_tresc = self.bezpieczny_tex(self._get_value(dane, "Treść tematu"))
        opis = self.bezpieczny_tex(self._get_value(dane, "Opis")).replace("\n", " ")
        jezyk = self.bezpieczny_tex(self._get_value(dane, "Język"))
        data_wydania = self.bezpieczny_tex(self._get_value(dane, "Data wydania"))
        zaklad = self.bezpieczny_tex(self._get_value(dane, "Zakład"))

        rodzaj_studiow = self.bezpieczny_tex(self._get_value(dane, "Rodzaj studiów"))
        if "niestacjonarne" in rodzaj_studiow:
            tex_rodzaj = r"\textbf{\sout{Stacjonarne} / Niestacjonarne (zaoczne) *}"
        else:
            tex_rodzaj = r"\textbf{Stacjonarne / \sout{Niestacjonarne (zaoczne)} *}"

        stopien_studiow = self.bezpieczny_tex(self._get_value(dane, "Stopień studiów"))
        if "magisterskie" in stopien_studiow or "ii" in stopien_studiow:
            tex_stopien = r"\textbf{\sout{I (inżynierskie)} / II (magisterskie) *}"
        else:
            tex_stopien = r"\textbf{I (inżynierskie) / \sout{II (magisterskie)} *}"

            
        
        status_tematu = self.bezpieczny_tex(self._get_value(dane, "Rodzaj studiów"))
        if "zmieniony" in status_tematu:
            tex_status = r"\noindent Temat pracy dyplomowej (\sout{nowy} / zmieniony)*:\\[0.1cm]"
        else:
            tex_status = r"\noindent Temat pracy dyplomowej (nowy / \sout{zmieniony})*:\\[0.1cm]"


        zadania_kody = []
        for i in range(1, 6):
            zad_kol = f"Zadanie {i}"
            if zad_kol in dane and pd.notna(dane[zad_kol]) and str(dane[zad_kol]).strip():
                zad_tekst = self.bezpieczny_tex(dane[zad_kol])
                zadania_kody.append(f"{i}. {zad_tekst}" + r" \dotfill\\[0.1cm]")

        if zadania_kody:
            zadania_kody[-1] = zadania_kody[-1].replace(r"\\[0.1cm]", r"\\[0.3cm]")
            
        # --- Budowa dokumentu (odwzorowanie oryginalnego układu) ---
        self.append(NoEscape(r"\setstretch{1.15}"))
        self.append(NoEscape(r"\pagestyle{empty}"))

        self.append(bold("Wydział Samochodów i Maszyn Roboczych"))

        self.append(NoEscape(r"\begin{center}"))
        self.append((NoEscape(r"\Large KARTA PRACY PRZEJŚCIOWEJ")))
        self.append(NoEscape(r"\end{center}"))
        #self.append(NewLine())

        # --- Dane podstawowe ---
        self.append(("Rodzaj studiów: "))
        self.append(f"{rodzaj_studiow}\n")

        self.append(("Stopień studiów: "))
        self.append(f"{stopien_studiow}\n")

        self.append(("Kierunek: "))
        self.append(f"{kierunek}\n")

        self.append(("Specjalność: "))
        self.append(f"{specjalnosc}\n\n")

        self.append(("Student: "))
        self.append(f"{dyplomant}  \n\n ")

        self.append(("Numer albumu: "))
        self.append(f"{album} \n\n")

        self.append(("Prowadzący: "))
        self.append(f"{prowadzacy} \n\n  ")

        self.append(("Zakład: "))
        self.append(f"{zaklad}\n")

        # --- Temat pracy ---
        status_tematu = "nowy"
        
        self.append((f"Temat pracy przejściowej ({status_tematu}):\n"))
        self.append(f"{temat_tresc}\n\n")

        # --- Opis ---
        self.append(("Zwięzły opis celu pracy:\n"))
        self.append(f"{opis}\n\n")

        # --- Daty i podpisy ---
        self.append(bold("Data wydania: "))
        self.append(f"{data_wydania}\n\n")

        self.append(bold("Data zwrotu**: "))
        self.append("09.06.2026\n\n")

        self.append(bold("Podpis studenta: "))
        self.append(".............................................................\n\n")

        self.append(bold("Podpis pracownika: "))
        self.append(".............................................................\n\n")

        # --- Informacja o zaliczeniu ---
        self.append(bold("INFORMACJA O ZALICZENIU\n\n"))

        self.append(bold("Data zaliczenia: "))
        self.append(".............................................................\n\n")

        self.append(bold("Ocena z pracy: "))
        self.append("5\n\n")

        self.append(bold("Podpis pracownika: "))
        self.append(".............................................................\n\n")

        # --- Przypisy ---
        self.append(
            NoEscape(
                r"** wpisać termin na tydzień przed końcem zajęć semestru."
            )
        )
        return self
    
    def sign_it(self, signature, role="promotor", offset=422, x_offset=152):
        import fitz
        import os

        pdf_bazowy = self.pdf_path
        pdf_wyjsciowy = f"{self.pdf_path}".replace('.pdf', '_signed.pdf')
        
        plik_wejsciowy = pdf_wyjsciowy if os.path.exists(pdf_wyjsciowy) else pdf_bazowy

        try:
            doc_pdf = fitz.open(plik_wejsciowy)
            strona = doc_pdf[0]

            if os.path.exists(signature):
                x1, y1 = x_offset, offset
                x2, y2 = x1 + 120, y1 + 40

                if role == "student":
                    # Tylko jeden podpis
                    prostokat_student = fitz.Rect(x1, y1, x2, y2)
                    strona.insert_image(prostokat_student, filename=signature)
                elif role == "promotor":
                    # Dwa podpisy (jeden niżej)
                    prostokat_promotor = fitz.Rect(x1, y1, x2, y2)
                    strona.insert_image(prostokat_promotor, filename=signature)

                    del_x, del_y = 10, 115
                    prostokat_promotor2 = fitz.Rect(x1 + del_x, y1 + del_y, x2 + del_x, y2 + del_y)
                    strona.insert_image(prostokat_promotor2, filename=signature)

                temp_pdf = pdf_wyjsciowy + ".tmp"
                doc_pdf.save(temp_pdf)
                doc_pdf.close()
                os.replace(temp_pdf, pdf_wyjsciowy)

        except Exception as e:
            print(f"Błąd przy wstawianiu podpisu: {e}")
    
import pandas as pd

#Tutaj wpisz swoje dane
dane_dict = {
    "Dyplomant": "Jan Kowalski",
    "Numer albumu": "123456",
    "Kierunek": "Mechanika i Budowa Maszyn",
    "Specjalność": "Pojazdy samochodowe",
    "Adres e-mail": "jan.kowalski@student.edu.pl",
    "Prowadzący": "dr inż. Nowak",
    "Konsultant": "mgr inż. Wiśniewski",
    "Treść tematu": "Analiza wytrzymałościowa ramy pojazdu",
    "Opis": "Celem pracy jest analiza wytrzymałościowa konstrukcji ramy pojazdu z wykorzystaniem MES.",
    "Język": "polski",
    "Data wydania": "01.10.2025",

    "Zadanie 1": "Przegląd literatury",
    "Zadanie 2": "Budowa modelu numerycznego",
    "Zadanie 3": "Przeprowadzenie symulacji",
    "Zadanie 4": "Analiza wyników",
    "Zadanie 5": "Opracowanie wniosków",
}

# część indeksowa (symuluje pd.Series z kolumnami bez nazw)
# indeks 1 -> rodzaj studiów
# indeks 2 -> stopień
# indeks 10 -> status tematu
dane_lista = [
    None,
    "stacjonarne",      # indeks 1
    "I stopnia",        # indeks 2
    None, None, None, None, None, None,
    None,
    "nowy"              # indeks 10
]


dane = pd.Series(dane_dict)
# for i, val in enumerate(dane_lista):
#     if val is not None:
#         dane.loc[i] = val

#Wywołanie
#doc = InterimProjectForm(
#   default_filepath="./output/karta_123456_Jan_Kowalski",
#   dane=dane_dict,
#   auto_build=True,
#   auto_generate=True  # wygeneruje .tex
#)
    pass
