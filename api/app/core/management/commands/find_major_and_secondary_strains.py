from django.core.management.base import BaseCommand
from core.models import Sample, SampleCount
from helpers.color_log import logger

class Command(BaseCommand):

    def handle(self, *args, **options):
        samples = Sample.objects.all()
        for sample in samples:
            sample.major_strain = None
            sample.secondary_strains.clear()
            sample.save()
            self.process_sample_counts(sample)


    def __get_dp10(self, dp):
        # dp string 1:0.00;2:0.00;3:0.00;4:0.00;5:0.00;6:0.00;7:0.00;8:0.00;9:0.00;10:0.00;11:0.00;12:0.00;13:0.00;14:0.00;15:0.00;16:0.00;17:0.00;18:0.00;19:0.00;20:0.00;
        dp_parts = dp.split(';')
        for part in dp_parts:
            if part.startswith('10:'):
                try:
                    return float(part.split(':')[1])
                except (ValueError, IndexError):
                    logger.warning(f"Error parsing DP10 value from part: {part}")
                    return 0.0
        logger.warning(f"DP10 value not found in DP string: {dp}")
        return 0.0

    def process_sample_counts(self, sample):
        sample_counts = SampleCount.objects.filter(sample=sample)
        filtered_viruses = []
        for sc in sample_counts:
            try:
                dp10 = self.__get_dp10(sc.DP) if sc.DP else 0.0
               # logger.info(f"Sample {sample.pseudonymized_id}, Substrain {sc.substrain.name}, Aligned {sc.aligned}, DP10 {dp10}, RPKM Proportions {sc.rpkm_proportions}")
                if sc.aligned and float(sc.aligned) > 10 and dp10 > 0.2:
                    filtered_viruses.append(sc)

            except (ValueError, TypeError) as e:
                logger.warning(
                    f"Error processing sample {sample.pseudonymized_id}, "
                    f"Substrain {getattr(sc.substrain, 'name', 'Unknown')}, "
                    f"DP {sc.DP}."
                    f" Exception: {e}"

                )
                continue

        secondary_strains = []
        major_strain_sc = None
        if len(filtered_viruses) == 1:
            major_strain_sc = filtered_viruses[0]
        elif len(filtered_viruses) > 1:
            filtered_viruses.sort(
                key=lambda x: (
                    self.__get_dp10(x.DP) if x.DP else 0.0,
                    x.rpkm_proportions if x.rpkm_proportions is not None else 0.0
                ),
                reverse=True
            )
            major_strain_sc = filtered_viruses[0]
            secondary_strains = filtered_viruses[1:]

        sample.major_strain = major_strain_sc.substrain if major_strain_sc else None
        sample.save()
        for strain in secondary_strains:
            sample.secondary_strains.add(strain.substrain)
        sample.save()
        if major_strain_sc or secondary_strains:
            logger.info(f"Sample {sample.pseudonymized_id} major strain: {sample.major_strain}, secondary strains: {[s.substrain.strain for s in secondary_strains]}")



# Per ogni sample nella tabella sample_counts cercare
# aligned_reads > 10
# DP10 > 0.2
# Questo filtro ti puo’ dare 0, 1 oppure >1 viruses
# Se 0 viruses, major strain e secondary_strain sono vuoti
# Se 1 virus, quel virus va in major strain
# If there is more than 1 virus with the same DP10, the major strain is the one with the highest rpkm_proportions, the others passing the filter are all secondary_strain
