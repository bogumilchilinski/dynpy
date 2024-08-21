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