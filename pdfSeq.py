from fpdf import FPDF

pdf = FPDF('P','in','Letter')
pdf.add_page()
pdf.set_margins(.5, .5)
pdf.set_font('Arial', 12)
pdf.cell(0, 0, 'Hello World!')
pdf.output('tuto1.pdf', 'F')