from pylatex import Document, Section, Subsection, Command, Itemize
from pylatex.utils import italic, NoEscape


def bool_to_str(bool_value):
    if bool_value:
        return 'Yes'
    else:
        return 'No'

def create_title(doc):
    doc.preamble.append(Command('title', 'Automatic Data Report'))
    doc.preamble.append(Command('author', 'R.H.'))
    doc.preamble.append(Command('date', NoEscape(r'\today')))
    doc.append(NoEscape(r'\maketitle'))
    
    
def create_information_secction(doc, dict_computation_options):
    with doc.create(Section('Computation resume')):
        doc.append('The electronic band strucutre was computed with Empirical Pseudopotential Method (EPM).\n')
        doc.append('The following options were used:')
        with doc.create(Itemize()) as itemize:
            for key, value in dict_computation_options.items():
                if type(value) is bool:
                    itemize.add_item(f"{key}: {bool_to_str(value)}")
                else:
                    itemize.add_item(f"{key}: {value}")

def fill_document(doc):
    """
    Fill the document with the different parts

    Parameters
    ----------
    doc : _type_
        _description_
    """
    create_title(doc)
    
    options = {'Method': 'EPM', "Number nearest neighbors": 6, "Number of bands": 6, "Number of kpoints": 2000, "path": "LKWGXWLGK",
               "Nonlocal Correction": True, "Spin Orbit Coupling": False}
    create_information_secction(doc, options)
        

def create_document(document_name):
    """Create a new document and fill it with content.

    :returns: the document
    :rtype: :class:`pylatex.document.Document` instance
    """
    doc = Document(document_name)
    return doc




if __name__ == '__main__':
    document_name = "data_report"
    doc = create_document(document_name)
    fill_document(doc)
    doc.generate_pdf('data_report', clean=False)