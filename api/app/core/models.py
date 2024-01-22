from django.db import models


# Create your models here.


class Plate(models.Model):
    barcode = models.TextField(unique=True)

    def __str__(self):
        return self.barcode


class Well(models.Model):
    plate = models.ForeignKey(Plate, on_delete=models.CASCADE)
    location = models.TextField(null=True)

    def __str__(self):
        return self.location


class Sample(models.Model):
    sample_number = models.TextField(unique=True)
    well = models.ForeignKey(Well, on_delete=models.CASCADE)
    pseudoanonymized_id = models.TextField(unique=True, null=True)
    plate = models.ForeignKey(Plate, on_delete=models.CASCADE, null=True)

    def __str__(self):
        return self.pseudoanonymized_id


class FileType(models.Model):
    postfix = models.TextField(unique=True)

    def __str__(self):
        return self.postfix


class File(models.Model):
    related_name = "files"
    path = models.TextField(unique=True)
    checksum = models.TextField(unique=True)
    type = models.ForeignKey(FileType, on_delete=models.CASCADE, null=True)
    sample = models.ForeignKey(Sample, on_delete=models.CASCADE, null=True)
    plate = models.ForeignKey(Plate, on_delete=models.CASCADE, null=True)

    def __str__(self):
        return self.path


class SampleFile(models.Model):
    sample = models.ForeignKey(Sample, on_delete=models.CASCADE)
    file_type = models.CharField(max_length=20)
    file = models.FileField(upload_to="uploads/")

    def __str__(self):
        return self.file.name


class Strain(models.Model):
    name = models.TextField()

    def __str__(self):
        return self.name


class Panel(models.Model):
    name = models.TextField()
    strain = models.ForeignKey(Strain, on_delete=models.CASCADE, null=True)

    def __str__(self):
        return self.name


class Substrain(models.Model):
    name = models.TextField()
    taxon_id = models.IntegerField(null=True)
    scientific_name = models.TextField(null=True)
    strain = models.ForeignKey(Strain, on_delete=models.CASCADE, null=True)

    def __str__(self):
        return self.name


class SampleCount(models.Model):
    plate = models.ForeignKey(Plate, on_delete=models.CASCADE, null=True)
    sample = models.ForeignKey(Sample, on_delete=models.CASCADE, null=True)
    substrain = models.ForeignKey(
        Substrain, on_delete=models.CASCADE, null=True
    )  # null=True because we first create the object by import and populate it later
    aligned = models.IntegerField(null=True)
    length = models.IntegerField(null=True)
    rpkm = models.FloatField(null=True)
    rpkm_proportions = models.FloatField(null=True)
    normcounts = models.FloatField(null=True)
    outlier = models.BooleanField(null=True)
    DP_threshold = models.FloatField(null=True, blank=True)
    DP = models.FloatField(null=True, blank=True)
    DP_status = models.CharField(max_length=20, null=True, blank=True)
    readnum_status = models.TextField(null=True, blank=True)
    readnum_threshold = models.FloatField(null=True, blank=True)
    percentile_threshold = models.TextField(null=True, blank=True)

    def __str__(self):
        return self.plate.barcode + " " + self.substrain.name


class Metadata(models.Model):
    plate = models.ForeignKey(Plate, on_delete=models.CASCADE, null=True)
    sample = models.ForeignKey(Sample, on_delete=models.CASCADE, null=True)
    well = models.ForeignKey(Well, on_delete=models.CASCADE, null=True)
    prescriber = models.CharField(max_length=20)
    order_date = models.DateField()
    ent_date = models.DateField()
    treatment_type = models.CharField(max_length=20, null=True)
    data = models.JSONField(null=True, default=list)
