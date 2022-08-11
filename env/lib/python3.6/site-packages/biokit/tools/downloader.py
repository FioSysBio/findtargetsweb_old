


def download_fasta(identifier, tag=None, service=None):
    """Retrieve FASTA from ENA, or UniProt, or EUtils

    :param identifier: valid identifier to retrieve FASTA (from ENA)
    :param tag. filename without extension (.fa)
    :param service: one of ENA, UniProt, EUtils to force the usage of one
        service (instead of 3). This is fasta since it avoids trying a service
        that may fail. See details below

    If the identifier is not found in ENA, we assume it is a protein sequence
    from uniprot. It it fails, we assume that the identifier is a NCBI
    identifier so we try EUtils service. If you know that the service is EUtils, 
    set service=EUtils which will skip the ENA and UniProt trials (that will
    fail).
    """
    from bioservices import ENA, UniProt, EUtils

    if service is None or service == "ENA":
        ena = ENA()
        data = ena.get_data(identifier, 'fasta')
        if isinstance(data, int) and data == 400 or "not supported" in data.decode():
            data = None

    # if not ENA, let us try UniProt
    if data is None and (service is None or service == "UniProt"):
        u = UniProt()
        data = u.get_fasta(identifier)
        if isinstance(data, int) and data == 400:
            data = None

    if data is None and (service is None or service == "EUtils"):
        ss = EUtils()
        data = ss.EFetch("nuccore", "808039399", rettype="fasta").decode()
        if isinstance(data, int) and data == 400:
            data = None


    if tag is None:
        return data
    else:
        with open("%s.fa" % tag, "w") as fout:
            fout.write(data)


def download_genbank(identifier, tag=None, db="nuccore", rettype="gbwithparts"):
    """

    :param identifier: valid identifier to retrieve from NCBI
    :param tag. filename without extension (.fa)

    ::

        from biokit.tools.downloader import download_genbank
        download_genbank('JB409847')

    """
    from bioservices import EUtils
    eu = EUtils()
    data = eu.EFetch(db=db, id=identifier, rettype=rettype, retmode="text")
    if isinstance(data, int) and data == 400:
        raise ValueError("%s not found on NCBI")

    if tag is None:
        return data.decode()
    else:
        with open("%s.gbk" %  tag, "w") as fout:
            fout.write(data.decode())
