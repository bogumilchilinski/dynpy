import re
from collections import Counter

class LatexAnalyzerV3:
    def __init__(self, latex_file, output_file=None):
        self.latex_file = latex_file
        self.output_file = output_file
        self.latex_content = None
        self.load_file()
    
    def load_file(self):
        try:
            with open(self.latex_file, 'r', encoding='utf-8') as file:
                self.latex_content = file.read()
        except FileNotFoundError:
            print(f"Error: {self.latex_file} not found.")

    def extract_text_from_commands(self, text):
        patterns = [
            r'\\section\*?\{([^}]*)\}',
            r'\\subsection\*?\{([^}]*)\}',
            r'\\subsubsection\*?\{([^}]*)\}',
            r'\\emph\{([^}]*)\}',
            r'\\caption\{([^}]*)\}'
        ]
        extracted_text = ''
        for pattern in patterns:
            matches = re.findall(pattern, text, re.IGNORECASE)
            extracted_text += ' '.join(matches) + ' '
        return extracted_text
    
    def strip_special_characters_and_commands(self, text):
        text_from_commands = self.extract_text_from_commands(text)
        cleaned_text = re.sub(r'\\[a-zA-Z]+\[[^\]]*\]\{[^}]*\}', '', text)  # Remove commands like \command[...] {...}
        cleaned_text = re.sub(r'\\[a-zA-Z]+\{[^}]*\}', '', cleaned_text)  # Remove commands like \command{...}
        cleaned_text = re.sub(r'\\[a-zA-Z]+', '', cleaned_text)  # Remove commands like \command
        cleaned_text = re.sub(r'[%$\\{}&_^#~]', '', cleaned_text)  # Remove special characters
        return text_from_commands + cleaned_text

    def count_characters(self):
        if self.latex_content:
            cleaned_text = self.strip_special_characters_and_commands(self.latex_content.lower())
            return len(cleaned_text)
        return 0

    def count_words(self, text):
        words = re.findall(r'\b\w+\b', text)
        return len(words)

    def top_words(self, n=10):
        if self.latex_content:
            cleaned_text = self.strip_special_characters_and_commands(self.latex_content.lower())
            words = re.findall(r'\b\w+\b', cleaned_text)
            word_counts = Counter(words)
            return word_counts.most_common(n)
        return []

    def count_words_and_characters_in_sections(self):
        if not self.latex_content:
            return {}

        content = self.latex_content
        sections = re.split(r'(\\section\*?\{[^}]*\})', content)
        counts = {}

        for i in range(1, len(sections), 2):
            section_title = sections[i]
            section_content = sections[i + 1]
            subsections = re.split(r'(\\subsection\*?\{[^}]*\})', section_content)

            section_words = 0
            section_characters = 0

            for j in range(1, len(subsections), 2):
                subsection_title = subsections[j]
                subsection_content = subsections[j + 1]
                subsubsections = re.split(r'(\\subsubsection\*?\{[^}]*\})', subsection_content)

                subsection_words = 0
                subsection_characters = 0

                for k in range(1, len(subsubsections), 2):
                    subsubsection_title = subsubsections[k]
                    subsubsection_content = subsubsections[k + 1]

                    cleaned_subsubsection_content = self.strip_special_characters_and_commands(subsubsection_content)
                    subsubsection_word_count = self.count_words(cleaned_subsubsection_content)
                    subsubsection_character_count = len(cleaned_subsubsection_content)

                    counts[subsubsection_title] = {
                        'words': subsubsection_word_count,
                        'characters': subsubsection_character_count
                    }

                    subsection_words += subsubsection_word_count
                    subsection_characters += subsubsection_character_count

                cleaned_subsection_content = self.strip_special_characters_and_commands(subsection_content)
                subsection_word_count = self.count_words(cleaned_subsection_content)
                subsection_character_count = len(cleaned_subsection_content)

                counts[subsection_title] = {
                    'words': subsection_word_count,
                    'characters': subsection_character_count
                }

                section_words += subsection_word_count
                section_characters += subsection_character_count

            cleaned_section_content = self.strip_special_characters_and_commands(section_content)
            section_word_count = self.count_words(cleaned_section_content)
            section_character_count = len(cleaned_section_content)

            counts[section_title] = {
                'words': section_word_count,
                'characters': section_character_count
            }

        return counts

class DynSysChecker:
    def __init__(self, system, debug= False, verbosity = False):
        self.system = system
        self.runnable = False
        self.debug = debug
        self.verbosity = verbosity 
        
    def check_init(self):

        if not hasattr(self.system, '__init__'):
            return False

        try:
            _ = self.system()
        except TypeError as e:
            print(e)
            return False
        self.runnable = True
        return True
    
    def check_components(self):

        if not hasattr(self.system, 'components'):
            return False
        
        if self.runnable:
            try:
                comp = self.system().components
            except Exception as e:
                if self.debug:
                    print(e)
                return False
        else:
            try:
                comp = self.system.components
            except Exception as e:
                print(e)
                return False
        if self.verbosity is True:
            return comp
        else:
            return True
    
    def check_default_data(self):

        if not hasattr(self.system, 'get_default_data'):
            return False
        
        if self.runnable:
            try:
                def_data = self.system().get_default_data()
            except Exception as e:
                if self.debug:
                    print(e)
                return False
        else:
            try:
                def_data = self.system.get_default_data()
            except Exception as e:
                return False
        if self.verbosity is True:
            return def_data
        else:
            return True
        
    def check_numerical_data(self):

        if not hasattr(self.system, 'get_numerical_data'):
            return False
        
        if self.runnable:
            try:
                data_num = self.system().get_numerical_data()
            except Exception as e:
                if self.debug:
                    print(e)
                return False
        else:
            try:
                data_num = self.system.get_numerical_data()
            except Exception:
                return False
        if self.verbosity is True:
            return data_num
        else:
            return True
        
    def check_symbols_description(self):

        if not hasattr(self.system, 'symbols_description'):
            return False
        
        if self.runnable:
            try:
                sym_desc = self.system().symbols_description()
            except Exception as e:
                if self.debug:
                    print(e)
                return False
        else:
            try:
                sym_desc = self.system.symbols_description()
            except Exception:
                return False
        if self.verbosity is True:
            return sym_desc
        else:
            return True
    
    def check_unit_dict(self):

        if not hasattr(self.system, 'unit_dict'):
            return False
        
        if self.runnable:
            try:
                unit_dict = self.system().unit_dict()
            except Exception as e:
                if self.debug:
                    print(e)
                return False
        else:
            try:
                unit_dict = self.system.unit_dict()
            except Exception:
                return False
        if self.verbosity is True:
            return unit_dict
        else:
            return True
    
    def check_equation(self):
        
        if self.runnable:
            try:
                equ = self.system()._eoms[0]
            except Exception as e:
                if self.debug:
                    print(e)
                return False
        else:
            try:
                equ = self.system._eoms[0]
            except Exception:
                return False
        if self.verbosity:
            return equ
        else:
            return True
        
    def check_picture(self):
        
        if self.runnable:
            try:
                self.system()._as_picture();
                pic = self.system().preview(example = 'detail_real_name');
            except Exception as e:
                if self.debug:
                    print(e)
                return False
        else:
            try:
                self.system._as_picture();
                pic = self.system.preview(example = 'detail_real_name');
            except Exception:
                return False
        if self.verbosity:
            return pic
        else:
            return True
    
class DynsysCheckerTable:
    def __init__(self, system, verbosity = False):
        import pandas as pd
        initial_data = {'System Name': [], 
                        'Init': [],
                        'Components': [],
                        'Default Data': [],
                        'Numerical Parameters': [],
                        'Symbols Description': [],
                        'Units': [],
                        'Equation': [],   
                        'Picture': [],
                        'Result':[]}
        
        self.df = pd.DataFrame(initial_data)
        pd.set_option("display.max_colwidth", None)
        pd.set_option('future.no_silent_downcasting', True)
        
        self.verbosity = verbosity
        
        if isinstance(system, __builtins__.__class__):
            self.list_from_Modulestructure(system)
        elif isinstance(system, str):
            self.list_from_Modulestructure(system)
        elif isinstance(system, list):
            self.table_loop(system)
        else:
            self.system = system
            self.table_append(system)
        
        
    def list_from_Modulestructure(self, path):
        from dynpy.utilities.creators import ModuleStructure
        lst = ModuleStructure(path).get_classes()

        import importlib  
        for element in lst:
            m = importlib.import_module(element[0])
            tmp = getattr(m, element[1])
            self.table_append(tmp()) 
    
    def table_append(self, sys):
            
            import pandas as pd
            
            try:
                self.df
            except:
                self.init_table()
            
            param = DynSysChecker(sys, verbosity=self.verbosity)
            
            new_row = {'System Name': sys.__class__.__name__, 
                        'Init': param.check_init(),
                        'Components': param.check_components(),
                        'Default Data': param.check_default_data(),
                        'Numerical Parameters': param.check_numerical_data(),
                        'Symbols Description': param.check_symbols_description(),
                        'Units': param.check_unit_dict(),
                        'Equation': param.check_equation(),   
                        'Picture': param.check_picture(),
                        'Result': []}
            if False in new_row.values():
                res = False
            else:
                res = True
            self.df = pd.concat([self.df, pd.DataFrame([new_row])], ignore_index=True)
            self.df.loc[self.df.index[-1], 'Result'] = res
            
   
    def table_loop(self, sys):
        for element in sys:
            self.table_append(element)
    
    def get_table(self):
        if self.verbosity is False:
            return self.df.replace({1.0: True, 0.0: False})
        return self.df