import os


def copy_pdf(source, destination):

    copy_comand = f" cp {source} {destination}"
    os.popen(copy_comand)


#     os.remove(f'{source}')
