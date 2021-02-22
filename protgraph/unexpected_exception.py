class UnexpectedException(Exception):

    def __init__(self, accession=None, aminoacid=None, position=None, message=None, additional_info=None):
        self.accession = accession
        self.aminoacid = aminoacid
        self.position = position
        self.message = message
        self.additional_info = additional_info
        super().__init__(message)

    def __str__(self):
        return "An unexpected case happend, which was not covered. Please report this issue by the maintainer.\n" \
               "Please provide the following Information:\n" \
               "Accession: {}, Aminoacid(s): {}, Position: {}\n" \
               "Additional Context: {}\n" \
               "Message: {}\n" \
               .format(self.accession, self.aminoacid, self.position, self.message, self.additional_info)
