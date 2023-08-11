import shortuuid


def generate_id(length=6):
    return shortuuid.ShortUUID().random(length=length)
